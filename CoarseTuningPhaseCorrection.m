function phi = CoarseTuningPhaseCorrection(MRS_struct, data, lipid_flag)

% Based on: Bao et al., 2013. A robust automatic phase correction method
% for signal dense spectra. J. Magn. Reson. 234, 82-89.
% doi:10.1016/j.jmr.2013.06.012

% ---------------- (1) Sub-experiment indexing ----------------

if ~MRS_struct.p.HERMES && length(MRS_struct.p.target) == 1 && any(strcmp(MRS_struct.p.target,{'GABAGlx','Lac','EtOH'}))
%     if ~lipid_flag % use NAA
%         ind = MRS_struct.fids.ON_OFF == 0;
%     else % use water
        ind = 1:size(data,2);
%     end
elseif ~MRS_struct.p.HERMES && length(MRS_struct.p.target) == 1 && strcmp(MRS_struct.p.target,'GSH')
    if ~lipid_flag % use NAA
        ind = 1:size(data,2);
    else
        % throw an error for now
    end
elseif MRS_struct.p.HERMES && ~MRS_struct.p.HERCULES && length(MRS_struct.p.target) == 2
    if ~lipid_flag % use NAA
        switch MRS_struct.p.vendor
            case 'GE'
                ind = [3:4:size(data,2) 4:4:size(data,2)];
            case {'Philips','Philips_data','Philips_raw'}
                ind = [1:4:size(data,2) 4:4:size(data,2)];
            case {'Siemens_twix','Siemens_rda','Siemens_dicom'}
                ind = [2:4:size(data,2) 3:4:size(data,2)];
        end
    else % use water
        switch MRS_struct.p.vendor
            case 'GE'
                ind = [2:4:size(data,2) 4:4:size(data,2)];
            case {'Philips','Philips_data','Philips_raw'}
                ind = [2:4:size(data,2) 4:4:size(data,2)];
            case {'Siemens_twix','Siemens_rda','Siemens_dicom'}
                ind = [2:4:size(data,2) 3:4:size(data,2)];
        end
    end
elseif MRS_struct.p.HERMES && ~MRS_struct.p.HERCULES && length(MRS_struct.p.target) == 3
    if ~lipid_flag % use NAA
        ind = [2:4:size(data,2) 3:4:size(data,2)];
    else % use water
        % throw an error for now
    end
elseif MRS_struct.p.HERMES && MRS_struct.p.HERCULES
    if ~lipid_flag % use NAA
        switch MRS_struct.p.vendor
            case {'Philips','Philips_data','Philips_raw'}
                % throw an error for now
            case {'Siemens_twix','Siemens_rda','Siemens_dicom'}
                ind = [3:4:size(data,2) 4:4:size(data,2)];
        end
    else % use water
        switch MRS_struct.p.vendor
            case {'Philips','Philips_data','Philips_raw'}
                % throw an error for now
            case {'Siemens_twix','Siemens_rda','Siemens_dicom'}
                ind = [2:4:size(data,2) 4:4:size(data,2)];
        end
    end
end

% ---------------- (2) Select which peak to phase ----------------

data = mean(data(:,ind),2);

ii = MRS_struct.ii;
freqRange = MRS_struct.p.sw(ii)/MRS_struct.p.LarmorFreq(ii);
freq = (size(data,1) + 1 - (1:size(data,1))) / size(data,1) * freqRange + 4.68 - freqRange/2;

k = 0.17;

if ~MRS_struct.p.HERMES && length(MRS_struct.p.target) == 1 && any(strcmp(MRS_struct.p.target,{'GABAGlx','Lac','EtOH'}))
    freqLim = freq <= 4.68+k & freq >= 4.68-k; % water
else
    if ~lipid_flag
        freqLim = freq <= 2.01+k & freq >= 2.01-k; % NAA
    else
        freqLim = freq <= 4.68+k & freq >= 4.68-k; % water
    end
end
[~,i] = max(abs(data(freqLim)));
% Shift freq range so max signal is centered
freq2 = freq(freqLim);
maxFreq = freq2(i);
freqLim = freq <= maxFreq+k & freq >= maxFreq-k;

in = data(freqLim);

% ---------------- (3) Phase correction ----------------

fminopts = optimset(@fminsearch);
fminopts = optimset(fminopts,'Display','off','TolFun',1e-6,'TolX',1e-6);
fun = @(x) objFunc(in,x);
phi = fminsearch(fun,0,fminopts);

if ishandle(45)
    close(45);
end

end


function out = objFunc(in,x)

phi = x(1);
in = real(in * exp(1i*phi));

% figure(45);
% plot(in);
% drawnow;

Start = sum(in(1:50));
End = sum(in(end-49:end));
out = abs(Start - End);

end



