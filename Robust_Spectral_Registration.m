function [AllFramesFTrealign, MRS_struct] = Robust_Spectral_Registration(MRS_struct)

% Robust spectral registration (MM: May 2019)

% Looping parameters
if MRS_struct.p.HERMES || MRS_struct.p.HERCULES % run registration four times - once for each HERMES/HERCULES sub-experiment
    SpecRegLoop = 3;
    SubspecToAlign = repmat([3 2 1 0], [1 size(MRS_struct.fids.data,2)/4]);
else % run registration twice for MEGA-PRESS acquisitions
    SpecRegLoop = 1;
    SubspecToAlign = MRS_struct.fids.ON_OFF;
end

ii = MRS_struct.ii;

% Pre-allocate memory
MRS_struct.out.SpecReg.freq(ii,:) = zeros(1,size(MRS_struct.fids.data,2));
MRS_struct.out.SpecReg.phase(ii,:) = zeros(1,size(MRS_struct.fids.data,2));
MRS_struct.out.SpecReg.MSE(ii,:) = zeros(1,size(MRS_struct.fids.data,2));
MRS_struct.fids.data_align = complex(zeros(size(MRS_struct.fids.data)));
DataToAlign = complex(zeros(size(MRS_struct.fids.data)));
count = 0;

% Inputs
time = (0:(MRS_struct.p.npoints(ii)-1))'/MRS_struct.p.sw(ii);
input.dwelltime = 1/MRS_struct.p.sw(ii);

% Optimization options
lsqnonlinopts = optimoptions(@lsqnonlin);
lsqnonlinopts = optimoptions(lsqnonlinopts,'Display','off','FunctionTolerance',1e-6,'OptimalityTolerance',1e-6,'StepTolerance',1e-6,'Algorithm','levenberg-marquardt');

% Automatic lipid/unstable residual water removal
freqRange = MRS_struct.p.sw(ii)/MRS_struct.p.LarmorFreq(ii);
freq = (MRS_struct.p.npoints(ii) + 1 - (1:MRS_struct.p.npoints(ii))) / MRS_struct.p.npoints(ii) * freqRange + 4.68 - freqRange/2;
lipidLim = freq <= 1.85 & freq >= 0;
noiseLim = freq <= 11 & freq >= 10;

S = mean(real(fftshift(fft(MRS_struct.fids.data,[],1),1)),2);
r = std(S(lipidLim)) / std(S(noiseLim));
q = std(MRS_struct.spec.F0freq(ii,:));

% if any(strcmp(MRS_struct.p.vendor,{'Siemens_twix','Siemens_rda','Siemens_dicom'}))
%     q_threshold = 0.04;
% else
%     q_threshold = 0.04;
% end

lipid_flag = 0;
unstable_H20_flag = 0;

if r > 40 || q > 0.04
    
    spec = fftshift(fft(MRS_struct.fids.data,[],1),1);
    
    if r > 40
        lipid_flag = 1;
    end
    if q > 0.04
        unstable_H20_flag = 1;
    end
    
    reverseStr = '';
    for jj = 1:size(MRS_struct.fids.data,2)
        
        if lipid_flag && ~unstable_H20_flag
            msg = sprintf('\nLipid contamination detected. Applying lipid filter to transient: %d\n', jj);
        elseif ~lipid_flag && unstable_H20_flag
            msg = sprintf('\nUnstable residual water detected. Applying residual water filter to transient: %d\n', jj);
        elseif lipid_flag && unstable_H20_flag
            msg = sprintf('\nLipid contamination and unstable residual water detected. Applying lipid and residual water filters to transient: %d\n', jj);
        end
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        
        DataToAlign(:,jj) = SignalFilter(spec(:,jj), lipid_flag, unstable_H20_flag, MRS_struct);
        
    end
    
    if ishandle(77)
        close(77);
    end
    
else
    
    DataToAlign = MRS_struct.fids.data;
    
end

% Spectral registration
while SpecRegLoop > -1
    
    % Use first n points of time-domain data, where n is the last point
    % where SNR > 2
    signal = abs(DataToAlign(:,SubspecToAlign == SpecRegLoop));
    noise = 2*std(signal(ceil(0.75*size(signal,1)):end,:));
    SNR = signal ./ repmat(noise, [size(DataToAlign,1) 1]);
    SNR = mean(SNR,2);
    n = find(SNR > 2);
    if isempty(n)
        tMax = 100;
    else
        tMax = n(end);
        if tMax < 100
            tMax = 100;
        end
    end
    
    % Flatten complex data for use in spectral registration
    clear flatdata;
    flatdata(:,1,:) = real(DataToAlign(1:tMax,SubspecToAlign == SpecRegLoop));
    flatdata(:,2,:) = imag(DataToAlign(1:tMax,SubspecToAlign == SpecRegLoop));
    
    % Determine optimal iteration order by calculating a similarity metric (mean squared error)
    D = zeros(size(flatdata,3));
    ind = find(SubspecToAlign == SpecRegLoop);
    for jj = 1:size(flatdata,3)
        for kk = 1:size(flatdata,3)
            tmp = sum((real(DataToAlign(1:tMax,ind(jj))) - real(DataToAlign(1:tMax,ind(kk)))).^2) / tMax;
            if tmp == 0
                D(jj,kk) = NaN;
            else
                D(jj,kk) = tmp;
            end
        end
    end
    d = nanmedian(D);
    [~,alignOrd] = sort(d);
    
    % Set initial reference transient based on similarity index
    flattarget = squeeze(flatdata(:,:,alignOrd(1)));
    target = flattarget(:);
    
    % Scalar to normalize transients (reduces optimization time)
    a = max(abs(target));
    
    % Pre-allocate memory
    if count == 0
        params = zeros(size(flatdata,3),2);
        MSE    = zeros(1,size(flatdata,3));
        w      = zeros(1,size(flatdata,3));
    end
    m = zeros(length(target),size(flatdata,3));
    
    % Starting values for optimization
    f0 = MRS_struct.spec.F0freq(ii,ind) * MRS_struct.p.LarmorFreq(ii);
    f0 = f0(alignOrd);
    f0 = f0 - f0(1);
    phi0 = zeros(size(f0));
    x0 = [f0(:) phi0(:)];
    
    % Determine frequency and phase offsets by spectral registration
    t = 0:input.dwelltime:(length(target)/2 - 1)*input.dwelltime;
    iter = 1;
    reverseStr = '';
    for jj = alignOrd
        
        msg = sprintf('\nRobust spectral registration - Iteration: %d', iter);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        
        transient = squeeze(flatdata(:,:,jj));
        input.data = transient(:)/a;
        
        fun = @(x) SpecReg(input, target/a, x);
        params(jj,:) = lsqnonlin(fun, x0(iter,:), [], [], lsqnonlinopts);
        
        f   = params(jj,1);
        phi = params(jj,2);
        m_c = complex(flatdata(:,1,jj), flatdata(:,2,jj));
        m_c = m_c .* exp(1i*pi*(t'*f*2+phi/180));
        m(:,jj) = [real(m_c); imag(m_c)];
        resid = target - m(:,jj);
        MSE(jj) = sum(resid.^2) / (length(resid) - 2);
        
        % Update reference
        w(jj) = 1 - sum((m(:,alignOrd(1)) - m(:,jj)).^2) / sum((m(:,alignOrd(1)) - mean(m(:,alignOrd(1)))).^2);
        if w(jj) < 0
            w(jj) = 0;
        end
        w(jj) = 0.5*w(jj);
        target = (1 - w(jj))*target + w(jj)*m(:,jj);
        
        iter = iter + 1;
        
    end
    
    count = count + 1;
    
    ind = find(SubspecToAlign == SpecRegLoop);
    
    MRS_struct.out.SpecReg.freq(ii,ind)  = params(:,1);
    MRS_struct.out.SpecReg.phase(ii,ind) = params(:,2);
    MRS_struct.out.SpecReg.MSE(ii,ind) = MSE;
    
    % Apply frequency and phase corrections
    for jj = 1:size(flatdata,3)
        MRS_struct.fids.data_align(:,ind(jj)) = MRS_struct.fids.data(:,ind(jj)) .* ...
            exp(1i*params(jj,1)*2*pi*time) * exp(1i*pi/180*params(jj,2));
    end
    
    if SpecRegLoop == 0
        
        % Align sub-spectra
        MRS_struct.fids.data_align = SubSpectralAlign(MRS_struct.fids.data_align, MRS_struct);
        
        % Line-broadening, zero-filling and FFT
        AllFramesFTrealign = MRS_struct.fids.data_align .* repmat((exp(-time*MRS_struct.p.LB*pi)), [1 size(MRS_struct.fids.data,2)]);
        AllFramesFTrealign = fftshift(fft(AllFramesFTrealign, MRS_struct.p.ZeroFillTo(ii),1),1);
        
        if ~MRS_struct.p.phantom
            % Zero-order phase correction
            phi = CoarseTuningPhaseCorrection(MRS_struct, AllFramesFTrealign, lipid_flag);
            AllFramesFTrealign = AllFramesFTrealign * exp(1i*phi);
            
            % Global frequency shift
            CrFreqRange = MRS_struct.spec.freq >= 2.925 & MRS_struct.spec.freq <= 3.125;
            [~,FrameMaxPos] = max(mean(real(AllFramesFTrealign(CrFreqRange,:)),2));
            freq = MRS_struct.spec.freq(CrFreqRange);
            CrFreqShift = freq(FrameMaxPos);
            CrFreqShift = CrFreqShift - 3.02;
            CrFreqShift_pts = round(CrFreqShift / abs(MRS_struct.spec.freq(1) - MRS_struct.spec.freq(2)));
            AllFramesFTrealign = circshift(AllFramesFTrealign, CrFreqShift_pts, 1);
        end
        
    end
    
    SpecRegLoop = SpecRegLoop - 1;
    
end

if ishandle(201)
    close(201);
end
if ishandle(200)
    close(200);
end
if ishandle(555)
    close(555);
end
fprintf('\n')

end



