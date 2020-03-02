function [AllFramesFTrealign, MRS_struct] = Robust_Spectral_Registration(MRS_struct)

% Robust spectral registration (MM: August 2019)

ii = MRS_struct.ii;

% Looping parameters
if MRS_struct.p.HERMES
    SpecRegLoop = 3;
    SubspecToAlign = repmat([3 2 1 0], [1 size(MRS_struct.fids.data,2)/4]);
else
    SpecRegLoop = 1;
    SubspecToAlign = MRS_struct.fids.ON_OFF;
end

% Pre-allocate memory
if MRS_struct.p.HERMES
    params = zeros(size(MRS_struct.fids.data,2)/4,2);
    MSE    = zeros(1,size(MRS_struct.fids.data,2)/4);
    w      = zeros(1,size(MRS_struct.fids.data,2)/4);
else
    params = zeros(size(MRS_struct.fids.data,2)/2,2);
    MSE    = zeros(1,size(MRS_struct.fids.data,2)/2);
    w      = zeros(1,size(MRS_struct.fids.data,2)/2);
end
MRS_struct.out.SpecReg.freq(ii,:)  = zeros(1,size(MRS_struct.fids.data,2));
MRS_struct.out.SpecReg.phase(ii,:) = zeros(1,size(MRS_struct.fids.data,2));
MRS_struct.out.SpecReg.MSE(ii,:)   = zeros(1,size(MRS_struct.fids.data,2));
MRS_struct.fids.data_align = complex(zeros(size(MRS_struct.fids.data)));
DataToAlign = complex(zeros(size(MRS_struct.fids.data)));

% Optimization options
lsqnonlinopts = optimoptions(@lsqnonlin);
lsqnonlinopts = optimoptions(lsqnonlinopts,'Algorithm','levenberg-marquardt','Display','off');

% Automatic unstable lipid/residual water removal
freqRange = MRS_struct.p.sw(ii)/MRS_struct.p.LarmorFreq(ii);
freq = (MRS_struct.p.npoints(ii) + 1 - (1:MRS_struct.p.npoints(ii))) / MRS_struct.p.npoints(ii) * freqRange + 4.68 - freqRange/2;
waterLim = freq <= 4.68 + 0.25 & freq >= 4.68 - 0.25;
lipidLim = freq <= 1.85 & freq >= 0;
noiseLim = freq <= 9 & freq >= 8;

S = mean(real(fftshift(fft(MRS_struct.fids.data,[],1),1)),2);
r = std(S(lipidLim)) / std(S(noiseLim));
r_threshold = 40;

spec = real(fftshift(fft(MRS_struct.fids.data,[],1),1));
if MRS_struct.p.HERMES
    ind = all(MRS_struct.fids.ON_OFF' == 0,2);
else
    switch MRS_struct.p.target{1}
        case {'GABAGlx','Lac','EtOH'}
            ind = 1:size(MRS_struct.fids.data,2);
        case 'GSH'
            ind = MRS_struct.fids.ON_OFF == 0;
    end
end
q = sum(abs(spec(waterLim,ind))) * abs(freq(1) - freq(2));
q = q / max(q);
q = sum(q < 0.5) / length(q);
q_threshold = 0.1;

lipid_flag = 0;
water_flag = 0;

if r > r_threshold || q > q_threshold
    if r > r_threshold
        lipid_flag = 1;
    end
    if q > q_threshold
        water_flag = 1;
    end
    
    spec = fftshift(fft(MRS_struct.fids.data,[],1),1);
    
    reverseStr = '';
    for jj = 1:size(MRS_struct.fids.data,2)
        if lipid_flag && ~water_flag
            msg = sprintf('\nUnstable lipid contamination detected. Applying lipid filter to transient: %d\n', jj);
        elseif ~lipid_flag && water_flag
            msg = sprintf('\nUnstable residual water detected. Applying residual water filter to transient: %d\n', jj);
        elseif lipid_flag && water_flag
            msg = sprintf('\nUnstable lipid contamination and residual water detected. Applying lipid and residual water filters to transient: %d\n', jj);
        end
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        DataToAlign(:,jj) = SignalFilter(spec(:,jj), lipid_flag, water_flag, MRS_struct);
    end
    if ishandle(77)
        close(77);
    end
else
    DataToAlign = MRS_struct.fids.data;
end

time = (0:(MRS_struct.p.npoints(ii)-1))'/MRS_struct.p.sw(ii);

% Spectral registration
while SpecRegLoop > -1
    
    % Use first n points of time-domain data, where n is the last point where abs(diff(mean(SNR))) > 0.5
    signal = abs(DataToAlign(:,SubspecToAlign == SpecRegLoop));
    noise = 2*std(signal(ceil(0.75*size(signal,1)):end,:));
    SNR = signal ./ repmat(noise, [size(DataToAlign,1) 1]);
    SNR = abs(diff(mean(SNR,2)));
    SNR = SNR(time <= 0.2); % use no more than 200 ms of data
    tMax = find(SNR > 0.5,1,'last');
    if isempty(tMax) || tMax < find(time <= 0.1,1,'last') % use at least 100 ms of data
        tMax = find(time <= 0.1,1,'last');
    end
    
    % Flatten complex data for use in spectral registration
    clear flatdata
    flatdata(:,1,:) = real(DataToAlign(1:tMax,SubspecToAlign == SpecRegLoop));
    flatdata(:,2,:) = imag(DataToAlign(1:tMax,SubspecToAlign == SpecRegLoop));
    
    % Determine optimal alignment order by calculating a similarity metric (mean squared error)
    if strcmp(MRS_struct.p.vendor,'Siemens_rda') % if .rda data, this subroutine doesn't apply
        alignOrd = 1;
    else
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
    end
    
    % Set initial reference transient based on similarity index
    target = squeeze(flatdata(:,:,alignOrd(1)));
    target = target(:);
    
    % Scalar to normalize transients (reduces optimization time)
    a = max(abs(target));
    
    % Pre-allocate memory
    m = zeros(length(target),size(flatdata,3));
    
    % Starting values for optimization
    f0 = MRS_struct.spec.F0freq2(ii,ind) * MRS_struct.p.LarmorFreq(ii);
    f0 = f0(alignOrd);
    f0 = f0 - f0(1);
    phi0 = zeros(size(f0));
    x0 = [f0(:) phi0(:)];
    
    % Determine frequency and phase offsets by spectral registration
    t = 0:(1/MRS_struct.p.sw(ii)):(length(target)/2-1)*(1/MRS_struct.p.sw(ii));
    iter = 1;
    reverseStr = '';
    for jj = alignOrd
        
        msg = sprintf('\nRobust spectral registration - Iteration: %d', iter);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        
        transient = squeeze(flatdata(:,:,jj));
        fun = @(x) SpecReg(transient(:)/a, target/a, t, x);
        params(jj,:) = lsqnonlin(fun, x0(iter,:), [], [], lsqnonlinopts);
        
        f   = params(jj,1);
        phi = params(jj,2);
        m_c = complex(flatdata(:,1,jj), flatdata(:,2,jj));
        m_c = m_c .* exp(1i*pi*(t'*f*2+phi/180));
        m(:,jj) = [real(m_c); imag(m_c)];
        resid = target - m(:,jj);
        MSE(jj) = sum(resid.^2) / (length(resid) - 2);
        
        % Update reference
        w(jj) = 0.5*corr(target, m(:,jj)).^2;
        target = (1 - w(jj))*target + w(jj)*m(:,jj);
        
        iter = iter + 1;
        
    end
    
    ind = find(SubspecToAlign == SpecRegLoop);    
    MRS_struct.out.SpecReg.freq(ii,ind)  = params(:,1);
    MRS_struct.out.SpecReg.phase(ii,ind) = params(:,2);
    MRS_struct.out.SpecReg.MSE(ii,ind) = MSE;
    
    % Apply frequency and phase corrections to raw data
    for jj = 1:size(flatdata,3)
        MRS_struct.fids.data_align(:,ind(jj)) = MRS_struct.fids.data(:,ind(jj)) .* ...
            exp(1i*params(jj,1)*2*pi*time) * exp(1i*pi/180*params(jj,2));
    end
    
    if SpecRegLoop == 0
        
        % Align subspectra
        MRS_struct.fids.data_align = SubSpectralAlign(MRS_struct.fids.data_align, water_flag, MRS_struct);
        
        % Line-broadening, zero-filling and FFT
        AllFramesFTrealign = MRS_struct.fids.data_align .* repmat((exp(-time*MRS_struct.p.LB*pi)), [1 size(MRS_struct.fids.data,2)]);
        AllFramesFTrealign = fftshift(fft(AllFramesFTrealign, MRS_struct.p.ZeroFillTo(ii), 1),1);
        
        if ~MRS_struct.p.phantom
            % Global frequency shift
            CrFreqRange = MRS_struct.spec.freq <= 3.02+0.15 & MRS_struct.spec.freq >= 3.02-0.15;
            [~,FrameMaxPos] = max(abs(mean(real(AllFramesFTrealign(CrFreqRange,:)),2)));
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



