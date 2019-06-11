function fids = SubSpectralAlign(fids, MRS_struct)

if MRS_struct.p.HERMES || MRS_struct.p.HERCULES
    data(:,1) = mean(fids(:,1:4:end),2);
    data(:,2) = mean(fids(:,2:4:end),2);
    data(:,3) = mean(fids(:,3:4:end),2);
    data(:,4) = mean(fids(:,4:4:end),2);
else
    data(:,1) = mean(fids(:,MRS_struct.fids.ON_OFF == 1),2);
    data(:,2) = mean(fids(:,MRS_struct.fids.ON_OFF == 0),2);
end

flatdata(:,1,:) = real(data);
flatdata(:,2,:) = imag(data);
data = abs(fftshift(fft(mean(data,2))));

ii = MRS_struct.ii;
freqRange = MRS_struct.p.sw(ii)/MRS_struct.p.LarmorFreq(ii);
freq = (size(data,1) + 1 - (1:size(data,1))) / size(data,1) * freqRange + 4.68 - freqRange/2;

if MRS_struct.p.HERMES || MRS_struct.p.HERCULES
    
    % Water
    freqLim(1,:) = freq <= 4.68+0.18 & freq >= 4.68-0.18;
    [~,i] = max(data(freqLim(1,:)));
    freq2 = freq(freqLim(1,:));
    maxFreq = freq2(i);
    freqLim(1,:) = freq <= maxFreq+0.18 & freq >= maxFreq-0.18;
    
    % NAA
    freqLim(2,:) = freq <= 2.01+0.13 & freq >= 2.01-0.13;
    [~,i] = max(data(freqLim(2,:)));
    freq2 = freq(freqLim(2,:));
    maxFreq = freq2(i);
    freqLim(2,:) = freq <= maxFreq+0.13 & freq >= maxFreq-0.13;
    
    % Cho
    freqLim(3,:) = freq <= 3.2+0.09 & freq >= 3.2-0.09;
    [~,i] = max(data(freqLim(3,:)));
    freq2 = freq(freqLim(3,:));
    maxFreq = freq2(i);
    freqLim(3,:) = freq <= maxFreq+0.06 & freq >= maxFreq-0.06;
    
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
    
else
    
    switch MRS_struct.p.target{1}
        case {'GABAGlx','Lac','EtOH'}
            % Water
            freqLim = freq <= freq <= 4.68+0.18 & freq >= 4.68-0.18;
            [~,i] = max(data(freqLim));
            freq2 = freq(freqLim);
            maxFreq = freq2(i);
            freqLim = freq <= maxFreq+0.18 & freq >= maxFreq-0.18;
        case 'GSH'
            % NAA
            freqLim = freq <= 2.01+0.13 & freq >= 2.01-0.13;
            [~,i] = max(data(freqLim));
            freq2 = freq(freqLim);
            maxFreq = freq2(i);
            freqLim = freq <= maxFreq+0.13 & freq >= maxFreq-0.13;
    end
    
end

dt = 1/MRS_struct.p.sw(ii);
t = 0:dt:(size(flatdata,1)-1)*dt;
lsqnonlinopts = optimoptions(@lsqnonlin);
lsqnonlinopts = optimoptions(lsqnonlinopts,'Display','off','FunctionTolerance',1e-6,'OptimalityTolerance',1e-6,'StepTolerance',1e-6);

if MRS_struct.p.HERMES || MRS_struct.p.HERCULES
    fun = @(x) objFunc(flatdata(:,:,[subSpecInd(2) subSpecInd(4)]), freqLim(1,:), t, x);
    param(1,:) = lsqnonlin(fun, [0 0], [], [], lsqnonlinopts);
    fun = @(x) objFunc(flatdata(:,:,[subSpecInd(1) subSpecInd(4)]), freqLim(2,:), t, x);
    param(2,:) = lsqnonlin(fun, [0 0], [], [], lsqnonlinopts);
    fun = @(x) objFunc(flatdata(:,:,[subSpecInd(3) subSpecInd(4)]), freqLim(3,:), t, x);
    param(3,:) = lsqnonlin(fun, [0 0], [], [], lsqnonlinopts);
else
    fun = @(x) objFunc(flatdata, freqLim, t, x);
    param = lsqnonlin(fun, [0 0], [], [], lsqnonlinopts);
end

if MRS_struct.p.HERMES || MRS_struct.p.HERCULES
    ind1 = subSpecInd(2):4:size(fids,2);
    ind2 = subSpecInd(1):4:size(fids,2);
    ind3 = subSpecInd(3):4:size(fids,2);
    for jj = 1:length(ind1)
        fids(:,ind1(jj)) = fids(:,ind1(jj)) .* exp(1i*param(1,1)*2*pi*t') * exp(1i*pi/180*param(1,2));
        fids(:,ind2(jj)) = fids(:,ind2(jj)) .* exp(1i*param(2,1)*2*pi*t') * exp(1i*pi/180*param(2,2));
        fids(:,ind3(jj)) = fids(:,ind3(jj)) .* exp(1i*param(3,1)*2*pi*t') * exp(1i*pi/180*param(3,2));
    end
else
    ind = find(MRS_struct.fids.ON_OFF == 1);
    for jj = 1:length(ind)
        fids(:,ind(jj)) = fids(:,ind(jj)) .* exp(1i*param(1)*2*pi*t') * exp(1i*pi/180*param(2));
    end
end

end


function out = objFunc(in, freqLim, t, x)

f   = x(1);
phi = x(2);

y = complex(in(:,1,1), in(:,2,1));
y = y .* exp(1i*pi*(t'*f*2+phi/180));

a = real(fftshift(fft(y)));
a = (a - min(a(freqLim))) / (max(a(freqLim)) - min(a(freqLim)));
b = real(fftshift(fft(complex(in(:,1,2), in(:,2,2)))));
b = (b - min(b(freqLim))) / (max(b(freqLim)) - min(b(freqLim)));

DIFF = a - b;

out = DIFF(freqLim);

end



