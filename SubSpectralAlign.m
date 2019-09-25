function fids = SubSpectralAlign(fids, water_flag, MRS_struct)

ii = MRS_struct.ii;
freqRange = MRS_struct.p.sw(ii)/MRS_struct.p.LarmorFreq(ii);
freq = (size(fids,1) + 1 - (1:size(fids,1))) / size(fids,1) * freqRange + 4.68 - freqRange/2;
time = (0:(MRS_struct.p.npoints(ii)-1))'/MRS_struct.p.sw(ii);
tMax = find(time <= 0.1,1,'last');

% Pre-allocate memory
if MRS_struct.p.HERMES
    n = 4;
else
    n = 2;
end
D = zeros(size(fids,2)/n);
w = cell(1,n);
data = complex(zeros(size(fids,1),n));
dataLim = [];

% Optimization options
lsqnonlinopts = optimoptions(@lsqnonlin);
lsqnonlinopts = optimoptions(lsqnonlinopts,'Algorithm','levenberg-marquardt','Display','off');

% Use weighted averaging to average subspectra
if strcmp(MRS_struct.p.vendor,'Siemens_rda') % if .rda data, use conventional averaging
    data = fids;
else
    for jj = 1:n
        if MRS_struct.p.HERMES
            ind = jj:n:size(fids,2);
        else
            ind = find(MRS_struct.fids.ON_OFF == abs(jj-2));
        end
        for kk = 1:size(fids,2)/n
            for ll = 1:size(fids,2)/n
                tmp = sum((real(fids(1:tMax,ind(kk))) - real(fids(1:tMax,ind(ll)))).^2) / 200;
                if tmp == 0
                    D(kk,ll) = NaN;
                else
                    D(kk,ll) = tmp;
                end
            end
        end
        d = nanmean(D);
        w{jj} = 1./d.^2;
        w{jj} = w{jj}/sum(w{jj});
        w{jj} = repmat(w{jj}, [size(fids,1) 1]);
        if water_flag
            dataLim = ceil(length(ind)/3);
            data(:,jj) = sum(w{jj}(:,1:dataLim) .* fids(:,ind(1:dataLim)),2);
        else
            data(:,jj) = sum(w{jj} .* fids(:,ind),2);
        end
    end
end
clear tmp
[~, data] = FlattenData(data);

% Set HERMES subexperiment indices (A, B, C, D) 
if MRS_struct.p.HERMES
    if ~MRS_struct.p.HERCULES
        if length(MRS_struct.p.target) == 2 && all(strcmp(MRS_struct.p.target,{'GABAGlx','GSH'}))
            switch MRS_struct.p.vendor
                case 'GE'
                    subSpecInd = [3 2 1 4];
                case {'Philips','Philips_data','Philips_raw'}
                    subSpecInd = [1 2 3 4];
                case {'Siemens_twix','Siemens_rda','Siemens_dicom'}
                    subSpecInd = [3 1 4 2];
            end
        elseif length(MRS_struct.p.target) == 3 && all(strcmp(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
            switch MRS_struct.p.vendor
                case {'Philips','Philips_data','Philips_raw'}
                    % throw an error for now
                case {'Siemens_twix','Siemens_rda','Siemens_dicom'}
                    subSpecInd = [3 1 4 2];
            end
        end
    else
        switch MRS_struct.p.vendor
            case {'Philips','Philips_data','Philips_raw'}
                subSpecInd = [1 4 3 2];
            case {'Siemens_twix','Siemens_rda','Siemens_dicom'}
                subSpecInd = [3 2 1 4];
        end
    end
end

% Phase-correct one subspectrum so all remaining subspectra are in the same phase
if MRS_struct.p.HERMES
    ind = subSpecInd(4);
else
    ind = 2;
end
fids = PhaseCorrection(data, fids, freq, ind, MRS_struct);

% Average subspectra again
if strcmp(MRS_struct.p.vendor,'Siemens_rda') % if .rda data, use conventional averaging
    data = fids;
else
    data = WeightedAveraging(fids, data, n, water_flag, dataLim, w, MRS_struct);
end
[flatdata, data] = FlattenData(data);

% Parameters for optimization
if MRS_struct.p.HERMES
    
    % Water
    freqLim(1,:) = freq <= 4.68+0.22 & freq >= 4.68-0.22;
    [~,i] = max(abs(data(freqLim(1,:), subSpecInd([2 4]))));
    freq2 = freq(freqLim(1,:));
    maxFreq = freq2(i);
    for jj = 1:2
        tmp(jj,:) = freq <= maxFreq(jj)+0.22 & freq >= maxFreq(jj)-0.22;
    end
    freqLim(1,:) = or(tmp(1,:), tmp(2,:));
    f0 = (maxFreq(1) - maxFreq(2)) * MRS_struct.p.LarmorFreq(ii);
    x0(1,:) = [f0 0];
    
    % NAA
    freqLim(2,:) = freq <= 2.01+0.13 & freq >= 2.01-0.13;
    [~,i] = max(abs(data(freqLim(2,:), subSpecInd([1 4]))));
    freq2 = freq(freqLim(2,:));
    maxFreq = freq2(i);
    for jj = 1:2
        tmp(jj,:) = freq <= maxFreq(jj)+0.18 & freq >= maxFreq(jj)-0.18;
    end
    freqLim(2,:) = or(tmp(1,:), tmp(2,:));
    f0 = (maxFreq(1) - maxFreq(2)) * MRS_struct.p.LarmorFreq(ii);
    x0(2,:) = [f0 0];
    
else
    
    switch MRS_struct.p.target{1}
        case {'GABAGlx','Lac','EtOH'}
            % Water
            freqLim = freq <= 4.68+0.22 & freq >= 4.68-0.22;
            [~,i] = max(abs(data(freqLim,:)));
            freq2 = freq(freqLim);
            maxFreq = freq2(i);
            for jj = 1:2
                freqLim(jj,:) = freq <= maxFreq(jj)+0.22 & freq >= maxFreq(jj)-0.22;
            end
%             % Cho
%             freqLim = freq <= 3.2+0.09 & freq >= 3.2-0.09;
%             [~,i] = max(abs(data(freqLim,:)));
%             freq2 = freq(freqLim);
%             maxFreq = freq2(i);
%             for jj = 1:2
%                 freqLim(jj,:) = freq <= maxFreq(jj)+0.08 & freq >= maxFreq(jj)-0.08;
%             end
        case 'GSH'
            % NAA
            freqLim = freq <= 2.01+0.13 & freq >= 2.01-0.13;
            [~,i] = max(abs(data(freqLim,:)));
            freq2 = freq(freqLim);
            maxFreq = freq2(i);
            for jj = 1:2
                freqLim(jj,:) = freq <= maxFreq(jj)+0.18 & freq >= maxFreq(jj)-0.18;
            end
    end
    freqLim = or(freqLim(1,:), freqLim(2,:));
    f0 = (maxFreq(1) - maxFreq(2)) * MRS_struct.p.LarmorFreq(ii);
    x0 = [f0 0];
    
end

t = 0:(1/MRS_struct.p.sw(ii)):(size(fids,1)-1)*(1/MRS_struct.p.sw(ii));

% Align subspectra to each other
if MRS_struct.p.HERMES
    
    % Water
    a = max(max(max(flatdata(:,:,subSpecInd([2 4])))));
    fun = @(x) objFunc(flatdata(:,:,subSpecInd([2 4]))./a, freqLim(1,:), t, x);
    param(1,:) = lsqnonlin(fun, x0(1,:), [], [], lsqnonlinopts);
    
    % NAA
    a = max(max(max(flatdata(:,:,subSpecInd([1 4])))));
    fun = @(x) objFunc(flatdata(:,:,subSpecInd([1 4]))./a, freqLim(2,:), t, x);
    param(2,:) = lsqnonlin(fun, x0(2,:), [], [], lsqnonlinopts);
    
    ind = subSpecInd(1):4:size(fids,2);
    for jj = 1:length(ind)
        fids(:,ind(jj)) = fids(:,ind(jj)) .* exp(1i*param(2,1)*2*pi*t') * exp(1i*pi/180*param(2,2));
    end
    for jj = 1:4
        ind = jj:4:size(fids,2);
        if water_flag
            data(:,jj) = sum(w{jj}(:,1:dataLim) .* fids(:,ind(1:dataLim)),2);
        else
            data(:,jj) = sum(w{jj} .* fids(:,ind),2);
        end
    end
    [flatdata, data] = FlattenData(data);
    
    % Cho
    freqLim(3,:) = freq <= 3.2+0.09 & freq >= 3.2-0.09;
    [~,i] = max(abs(data(freqLim(3,:), subSpecInd([3 1]))));
    freq2 = freq(freqLim(3,:));
    maxFreq = freq2(i);
    for jj = 1:2
        tmp(jj,:) = freq <= maxFreq(jj)+0.08 & freq >= maxFreq(jj)-0.08;
    end
    freqLim(3,:) = or(tmp(1,:), tmp(2,:));
    f0 = (maxFreq(1) - maxFreq(2)) * MRS_struct.p.LarmorFreq(ii);
    x0(3,:) = [f0 0];
    
    a = max(max(max(flatdata(:,:,subSpecInd([3 1])))));
    fun = @(x) objFunc(flatdata(:,:,subSpecInd([3 1]))./a, freqLim(3,:), t, x);
    param(3,:) = lsqnonlin(fun, x0(3,:), [], [], lsqnonlinopts);
    
else
    
    a = max(max(max(flatdata)));
    switch MRS_struct.p.target{1}
        case {'GABAGlx','Lac','EtOH'}
            fun = @(x) objFunc(flatdata./a, freqLim, t, x);
        case 'GSH'
            fun = @(x) objFunc(flatdata./a, freqLim, t, x);
    end
    param = lsqnonlin(fun, x0, [], [], lsqnonlinopts);
    
end

% Apply frequency and phase corrections to subspectra
if MRS_struct.p.HERMES
    ind1 = subSpecInd(2):4:size(fids,2);
    ind2 = subSpecInd(3):4:size(fids,2);
    for jj = 1:length(ind1)
        fids(:,ind1(jj)) = fids(:,ind1(jj)) .* exp(1i*param(1,1)*2*pi*t') * exp(1i*pi/180*param(1,2));
        fids(:,ind2(jj)) = fids(:,ind2(jj)) .* exp(1i*param(3,1)*2*pi*t') * exp(1i*pi/180*param(3,2));
    end
else
    ind = find(MRS_struct.fids.ON_OFF == 1);
    for jj = 1:length(ind)
        fids(:,ind(jj)) = fids(:,ind(jj)) .* exp(1i*param(1)*2*pi*t') * exp(1i*pi/180*param(2));
    end
end

if ishandle(44)
    close(44);
end
if ishandle(3)
    close(3);
end

end


function [flatdata, data] = FlattenData(data)

flatdata(:,1,:) = real(data);
flatdata(:,2,:) = imag(data);
data = real(fftshift(fft(data,[],1),1));

end


function fids = PhaseCorrection(data, fids, freq, ind, MRS_struct)

OFF = data(:,ind);
ii = MRS_struct.ii;

freqLim = freq <= 3.02+0.15 & freq >= 3.02-0.15;
[~,i] = max(abs(OFF(freqLim)));
freq2 = freq(freqLim);
maxFreq = freq2(i);
freqLim = freq <= maxFreq+0.58 & freq >= maxFreq-0.42;
OFF = OFF(freqLim);
Baseline = (OFF(1) + OFF(end))/2;
Width = 0.05;
Area = (max(OFF) - min(OFF)) * Width * 4;

x0 = [Area Width maxFreq 0 Baseline 0 1] .* [1 2*MRS_struct.p.LarmorFreq(ii) MRS_struct.p.LarmorFreq(ii) 180/pi 1 1 1];
ModelParamChoCr = FitChoCr(freq(freqLim), OFF, x0, MRS_struct.p.LarmorFreq(ii));

fids = fids * exp(1i*pi/180*ModelParamChoCr(4));

end


function data = WeightedAveraging(fids, data, n, water_flag, dataLim, w, MRS_struct)

for ii = 1:n
    if MRS_struct.p.HERMES
        ind = ii:4:size(fids,2);
    else
        ind = find(MRS_struct.fids.ON_OFF == abs(ii-2));
    end
    if water_flag
        data(:,ii) = sum(w{ii}(:,1:dataLim) .* fids(:,ind(1:dataLim)),2);
    else
        data(:,ii) = sum(w{ii} .* fids(:,ind),2);
    end
end

end


function out = objFunc(in, freqLim, t, x)

f   = x(1);
phi = x(2);

y = complex(in(:,1,1), in(:,2,1));
y = y .* exp(1i*pi*(t'*f*2+phi/180));

a = real(fftshift(fft(y)));
b = real(fftshift(fft(complex(in(:,1,2), in(:,2,2)))));

DIFF = a - b;
out = DIFF(freqLim);

end



