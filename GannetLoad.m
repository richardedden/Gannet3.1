function MRS_struct = GannetLoad(varargin)
% Gannet 3.1 GannetLoad
% Started by RAEE Nov. 5, 2012
% Updates by MGS, MM, GO 2016-2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Workflow summary
%   1. Pre-initialise
%   2. Determine data parameters from headers
%   3. Load data from files
%   4. Reconstruction of coil-sensitivity maps (PRIAM only)
%   5. Apply appropriate pre-processing
%   6. Output processed spectra
%   7. Build GannetLoad output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MRS_struct.version.Gannet = '3.1.1';
MRS_struct.version.load = '190612'; % set to date when final updates have been made

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   0. Check the file list for typos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metabfile = varargin{1};
missing = 0;
for filecheck = 1:length(metabfile)
    % If only water-suppressed data are provided, select Cr as reference.
    MRS_struct.p.Reference_compound = 'Cr';
    if ~exist(metabfile{filecheck},'file')
        disp(['The file ' metabfile{filecheck} ' (' num2str(filecheck) ')' ' is missing. Typo?']);
        missing = 1;
    end
end

if nargin < 3
    mode = 'batch';
    if nargin < 2
    else
        if iscell(varargin{2})
            % If water-unsuppressed data are provided, select H2O as reference.
            waterfile = varargin{2};
            MRS_struct.waterfile = waterfile;
            MRS_struct.p.Reference_compound = 'H2O';
            for filecheck = 1:length(waterfile)
                if ~exist(waterfile{filecheck},'file')
                    disp(['The file ' waterfile(filecheck) ' is missing. Typo?'])
                    missing = 1;
                end
            end
        else
            mode = varargin{2};
        end
    end
else
    mode = varargin{3};
    if ~strcmp(mode, 'batch') && ~strcmp(mode, 'join')
        error('The third input argument needs to be either ''batch'' or ''join''.')
    end
    % If water-unsuppressed data are provided, select H2O as reference.
    waterfile = varargin{2};
    MRS_struct.waterfile = waterfile;
    MRS_struct.p.Reference_compound = 'H2O';
    for filecheck = 1:length(waterfile)
        if ~exist(waterfile{filecheck},'file')
            disp(['The file ' waterfile(filecheck) ' is missing. Typo?'])
            missing = 1;
        end
    end
end

if missing
    error('Not all the files are there, so I give up.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1. Pre-initialise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MRS_struct.ii = 0;
MRS_struct.metabfile = metabfile;
MRS_struct = GannetPreInitialise(MRS_struct);
CheckTargets(MRS_struct); 

if MRS_struct.p.PRIAM
    vox = MRS_struct.p.Vox;
else
    vox = MRS_struct.p.Vox(1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2. Determine data parameters from header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Discern input data format
MRS_struct = GannetDiscernDatatype(metabfile{1}, MRS_struct);

% Determine number of provided water-suppressed files in the batch
MRS_struct.p.Reference_compound = 'Cr';
switch mode
    case 'batch'
        numscans = numel(metabfile);
        % For Siemens RDA, each acquisition has two RDA files, i.e. correct the
        % number:
        if strcmpi(MRS_struct.p.vendor,'Siemens_rda')
            numscans = numscans/2;
        end
        numfilesperscan = 1;
    case 'join'
        numscans = 1;
        numfilesperscan = numel(metabfile);
        % For Siemens RDA, each acquisition has two RDA files, i.e. correct the
        % number:
        if strcmpi(MRS_struct.p.vendor,'Siemens_rda')
            numfilesperscan = numfilesperscan/2;
        end
        fprintf('Running GannetLoad in ''join'' mode. FIDs from %i separate files...\n', numfilesperscan);
end

% Determine number of provided water-unsuppressed files in the batch
if exist('waterfile','var')
    MRS_struct.p.Reference_compound = 'H2O';
    numwaterscans = numel(waterfile);
    switch mode
        case 'batch'
            if numwaterscans ~= numscans
                error ('Number of water-unsuppressed files does not match number of water-suppressed files.');
            end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   3. Load data from files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MRS_struct.p.numscans = numscans;

for ii = 1:numscans % Loop over all files in the batch (from metabfile)
    
    MRS_struct.ii = ii;
    
    switch MRS_struct.p.vendor
        
        case 'GE'
            
            MRS_struct = GERead(MRS_struct, metabfile{ii});
            WaterData = MRS_struct.fids.data_water;
            MRS_struct.p.Reference_compound = 'H2O';
            MRS_struct.fids.data = MRS_struct.fids.data * MRS_struct.p.nrows(ii)/MRS_struct.p.Navg(ii);
            FullData = MRS_struct.fids.data;
            % Determine order of ON and OFF acquisitions
            MRS_struct = SpecifyOnOffOrder(MRS_struct);
            
        case 'Siemens_twix'
            
            if exist('waterfile','var')
                if strcmp(mode,'batch')
                    MRS_struct = SiemensTwixRead(MRS_struct, metabfile{ii}, waterfile{ii});
                else
                    % Load each input file and append the FIDs
                    MRS_struct = SiemensTwixRead(MRS_struct, metabfile{ii}, waterfile{ii});
                    for kk = 2:numfilesperscan
                        sub_MRS_struct = SiemensTwixRead(MRS_struct, metabfile{kk}, waterfile{ii});
                        MRS_struct.fids.data = [MRS_struct.fids.data sub_MRS_struct.fids.data];
                    end
                end
                % Correct the total number of averages
                MRS_struct.p.nrows = MRS_struct.p.nrows * numfilesperscan;
                MRS_struct.p.Navg = MRS_struct.p.Navg * numfilesperscan;
                WaterData = MRS_struct.fids.data_water;
            else
                if strcmp(mode,'batch')
                    MRS_struct = SiemensTwixRead(MRS_struct, metabfile{ii});
                else
                    % Load each input file and append the FIDs
                    MRS_struct = SiemensTwixRead(MRS_struct, metabfile{1});
                    for kk = 2:numfilesperscan
                        sub_MRS_struct = SiemensTwixRead(MRS_struct, metabfile{kk});
                        MRS_struct.fids.data = [MRS_struct.fids.data sub_MRS_struct.fids.data];
                    end
                end
                % Correct the total number of averages
                MRS_struct.p.nrows = MRS_struct.p.nrows * numfilesperscan;
                MRS_struct.p.Navg = MRS_struct.p.Navg * numfilesperscan;
            end
            % MM (160914): Need to set Water_Positive based on water signal
            if MRS_struct.p.Water_Positive == 0
                MRS_struct.fids.data = -MRS_struct.fids.data;
            end
            % Determine order of ON and OFF acquisitions
            MRS_struct = SpecifyOnOffOrder(MRS_struct);
            FullData = MRS_struct.fids.data;
            
        case 'Siemens_dicom' % GO 11/01/2016
            
            if exist('waterfile','var')
                MRS_struct = SiemensDICOMRead(MRS_struct,metabfile{ii},waterfile{ii}); % GO 02/05/2017
                WaterData = MRS_struct.fids.data_water;
            else
                MRS_struct = SiemensDICOMRead(MRS_struct,metabfile{ii}); % GO 11/01/2016
            end
            FullData = MRS_struct.fids.data;
            % Determine order of ON and OFF acquisitions
            MRS_struct = SpecifyOnOffOrder(MRS_struct);
            
        case 'dicom' % GO 11/30/2016
            
            if exist('waterfile','var')
                MRS_struct = DICOMRead(MRS_struct,metabfile{ii},waterfile{ii}); % GO 02/05/2017
                WaterData = MRS_struct.fids.data_water;
            else
                MRS_struct = DICOMRead(MRS_struct,metabfile{ii}); % GO 11/01/2016
            end
            FullData = MRS_struct.fids.data;
            
            % fill up fields required for downstream processing % GO 11/30/2016
            switch MRS_struct.p.ONOFForder
                % Not sure whether this is always the case, but the CMRR
                % sequence appears to go OFF-OFF-ON-ON in the DICOM
                % sorting?! Fixing this hard for now. GO 112017
                case 'onfirst'
%                     if strcmp(MRS_struct.p.seq,'""%CustomerSeq%\eja_svs_mpress""')
%                         MRS_struct.fids.ON_OFF=repmat([1 1 0 0],[1 MRS_struct.p.Navg(ii)/4]);
%                         MRS_struct.fids.ON_OFF=MRS_struct.fids.ON_OFF(:).';
%                     else
                        MRS_struct.fids.ON_OFF=repmat([1 0],[1 MRS_struct.p.Navg(ii)/2]);
                        MRS_struct.fids.ON_OFF=MRS_struct.fids.ON_OFF(:).';
%                    end
                case 'offfirst'
%                     if strcmp(MRS_struct.p.seq,'""%CustomerSeq%\eja_svs_mpress""')
%                         MRS_struct.fids.ON_OFF=repmat([0 0 1 1],[1 MRS_struct.p.Navg(ii)/4]);
%                         MRS_struct.fids.ON_OFF=MRS_struct.fids.ON_OFF(:).';
%                     else
                        MRS_struct.fids.ON_OFF=repmat([0 1],[1 MRS_struct.p.Navg(ii)/2]);
                        MRS_struct.fids.ON_OFF=MRS_struct.fids.ON_OFF(:).';
%                    end
            end
            
        case 'Siemens_rda'
            
            if exist('waterfile','var')
                if strcmp(mode,'batch')
                    MRS_struct = SiemensRead(MRS_struct, metabfile{ii*2}, metabfile{ii*2-1}, waterfile{ii});
                else
                    % Load each input file and append the FIDs
                    MRS_struct = SiemensRead(MRS_struct, metabfile{2}, metabfile{1}, waterfile{ii});
                    for kk = 2:numfilesperscan
                        sub_MRS_struct = SiemensRead(MRS_struct, metabfile{kk*2}, metabfile{kk*2-1}, waterfile{ii});
                        MRS_struct.fids.data = [MRS_struct.fids.data sub_MRS_struct.fids.data];
                    end
                end
                % Correct the total number of averages
                MRS_struct.p.Navg = MRS_struct.p.Navg * numfilesperscan;
                WaterData = MRS_struct.fids.data_water;
                MRS_struct.p.Nwateravg = 1;
            else
                if strcmp(mode,'batch')
                    MRS_struct = SiemensRead(MRS_struct, metabfile{ii*2}, metabfile{ii*2-1});
                else
                    % Load each input file and append the FIDs
                    MRS_struct = SiemensRead(MRS_struct, metabfile{2}, metabfile{1});
                    for kk = 2:numfilesperscan
                        sub_MRS_struct = SiemensRead(MRS_struct, metabfile{kk*2}, metabfile{kk*2-1});
                        MRS_struct.fids.data = [MRS_struct.fids.data sub_MRS_struct.fids.data];
                    end
                end
                % Correct the total number of averages
                MRS_struct.p.Navg = MRS_struct.p.Navg * numfilesperscan;
            end
            FullData = MRS_struct.fids.data;
            % Determine order of ON and OFF acquisitions
            MRS_struct = SpecifyOnOffOrder(MRS_struct);
            
        case 'Philips'
            
            if exist('waterfile','var')
                MRS_struct = PhilipsRead(MRS_struct, metabfile{ii}, waterfile{ii});
                WaterData = MRS_struct.fids.data_water;
            else
                MRS_struct = PhilipsRead(MRS_struct, metabfile{ii});
            end
            % Need to set Water_Positive based on water signal
            if MRS_struct.p.Water_Positive == 0
                MRS_struct.fids.data = -MRS_struct.fids.data;
            end
            FullData = MRS_struct.fids.data;
            MRS_struct = SpecifyOnOffOrder(MRS_struct);
            
        case 'Philips_data'
            
            % If a water reference scan is acquired, it is saved as a mix
            % in the DATA/LIST files. Later: add option to provide an additional
            % water reference file (i.e. short-TE). GO 03/02/2018
            if nargin > 1
                MRS_struct = PhilipsRead_data(MRS_struct, metabfile{ii}, waterfile{ii});
            else
                MRS_struct = PhilipsRead_data(MRS_struct, metabfile{ii});
            end
            if isfield(MRS_struct.fids, 'data_water')
                MRS_struct.p.Reference_compound = 'H2O';
                WaterData = MRS_struct.fids.data_water;
            else
                MRS_struct.p.Reference_compound = 'Cr';
            end
            MRS_struct = SpecifyOnOffOrder(MRS_struct); %For 3T and 7T -- 08212018 MGSaleh
            FullData = MRS_struct.fids.data;
        
        case 'Philips_raw' % GO 11/01/2016
            
            MRS_struct = PhilipsRawLoad(MRS_struct,metabfile{ii},3,0); % GO 11/02/2016 
            MRS_struct.fids.data=conj(squeeze(MRS_struct.multivoxel.allsignals(:,:,1,:)));
            if exist('waterfile','var')
                MRS_struct.p.Reference_compound = 'H2O';
                WaterData = MRS_struct.fids.data_water; % GO 11/03/2016
            end % GO 11/03/2016
            FullData = MRS_struct.fids.data;
            MRS_struct = SpecifyOnOffOrder(MRS_struct);
    
    end % end of vendor switch loop for data load
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   4. Reconstruction of coil-sensitivity maps
    %      (PRIAM only)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if a PRIAM dataset is processed, load the coil reference scan and
    % calculate the SENSE reconstruction matrix here
    if MRS_struct.p.PRIAM
        MRS_struct = senseRecon(MRS_struct);
        PRIAMData = zeros(length(MRS_struct.p.Vox),MRS_struct.p.Navg,MRS_struct.p.npoints);
        PRIAMWaterData = zeros(length(MRS_struct.p.Vox),MRS_struct.p.Nwateravg,MRS_struct.p.npoints);
        for kk = 1:MRS_struct.p.Navg
            PRIAMData(:,kk,:) = MRS_struct.p.SENSE.U * squeeze(FullData(:,kk,:));
            % Phase by multiplying with normalized complex conjugate of first point
            conj_norm = conj(PRIAMData(:,kk,1)) ./ abs(conj(PRIAMData(:,kk,1)));
            PRIAMData(:,kk,:) = PRIAMData(:,kk,:) .* repmat(conj_norm, [1 1 MRS_struct.p.npoints]);
        end
        for kk = 1:MRS_struct.p.Nwateravg
            PRIAMWaterData(:,kk,:) = MRS_struct.p.SENSE.U * squeeze(WaterData(:,kk,:));
            % Phase by multiplying with normalized complex conjugate of first point
            conj_norm = conj(PRIAMWaterData(:,kk,1)) ./ abs(conj(PRIAMWaterData(:,kk,1)));
            PRIAMWaterData(:,kk,:) = PRIAMWaterData(:,kk,:) .* repmat(conj_norm, [1 1 MRS_struct.p.npoints]);
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   5. Apply appropriate pre-processing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for kk = 1:length(vox) % loop over number of voxels % GO 03/26/2018
        
        % Select data from first voxel
        if MRS_struct.p.PRIAM
            FullData = squeeze(-PRIAMData(kk,:,:))';
            MRS_struct.fids.data = FullData;
            WaterData = squeeze(PRIAMWaterData(kk,:,:))';
            MRS_struct.fids.data_water = WaterData;
        end
        
        % Zero-fill to obtain nominal spectral resolution of 0.061 Hz/point
        MRS_struct.p.ZeroFillTo(ii) = round(32768/2000*MRS_struct.p.sw(ii));
        MRS_struct.p.zf = MRS_struct.p.ZeroFillTo(ii)/MRS_struct.p.npoints(ii);
        time = (1:1:size(FullData,1))/MRS_struct.p.sw(ii);
        
        % Finish processing water data
        if strcmpi(MRS_struct.p.Reference_compound,'H2O')
            
            if strcmpi(MRS_struct.p.vendor,'GE')
                ComWater = mean(WaterData,2);
            elseif strcmpi(MRS_struct.p.vendor,'Siemens_rda')
                ComWater = WaterData;
            elseif strcmpi(MRS_struct.p.vendor,'Siemens_twix')
                ComWater = WaterData;
            elseif (strcmpi(MRS_struct.p.vendor,'Siemens_dicom')) % GO 02/05/2017
                ComWater = mean(WaterData,2);
            elseif (strcmpi(MRS_struct.p.vendor,'dicom')) % GO 02/05/2017
                ComWater = mean(WaterData,2);
            elseif (strcmpi(MRS_struct.p.vendor,'Philips_raw')) % GO 02/05/2017
                ComWater = mean(WaterData(kk,:,:),2);
            elseif (strcmpi(MRS_struct.p.vendor,'Philips_data')) % GO 03/18/2018
                ComWater = mean(WaterData,2);
            else
                ComWater = WaterData.';
            end
        
            % Performing phase corrrection on the water-suppressed data
            % based on Klose (1990), MRM,14:26-30. The equation was
            % taken from Jiru (2008), EJR,67:202-217 -- MGSaleh 2016
            if MRS_struct.p.data_phase_correction
                if any(strcmpi(MRS_struct.p.vendor,{'Philips','Philips_data'}))
                    MRS_struct.fids.data = phase_correction_fids(MRS_struct.fids.data.', ComWater.');
                    MRS_struct.fids.data = MRS_struct.fids.data.';
                    FullData = MRS_struct.fids.data;
                else
                    MRS_struct.fids.data = phase_correction_fids(MRS_struct.fids.data.', ComWater);
                    MRS_struct.fids.data = MRS_struct.fids.data.';
                    FullData = MRS_struct.fids.data;
                end
            end
            
            % Performing phase corrrection on the unsuppressed water data
            if MRS_struct.p.water_phase_correction
                if any(strcmpi(MRS_struct.p.vendor,{'Philips','Philips_data'}))
                    ComWater = phase_correction_fids(ComWater.', ComWater.');
                    ComWater = ComWater.';
                else
                    ComWater = phase_correction_fids(ComWater, ComWater);
                end
            end
            
            % Line-broadening, zero-filling and FFT
            % GO (180514): Water data may have different bandwidth
            if isfield(MRS_struct.p,'sw_water')
                time_water = (1:1:size(ComWater,1))/MRS_struct.p.sw_water(ii);
            else
                time_water = (1:1:size(ComWater,1))/MRS_struct.p.sw(ii);
            end
            ComWater = ComWater .* exp(-time_water'*MRS_struct.p.LB*pi);
            MRS_struct.spec.(vox{kk}).water(ii,:) = fftshift(fft(ComWater,MRS_struct.p.ZeroFillTo(ii),1),1)';
            
        end % end of H2O reference loop
        
        % Line-broadening, zero-filling and FFT
        FullData = FullData .* repmat((exp(-time'*MRS_struct.p.LB*pi)), [1 size(FullData,2)]);
        MRS_struct.fids.FullData = FullData;
        AllFramesFT = fftshift(fft(FullData,MRS_struct.p.ZeroFillTo(ii),1),1);
        
        % Work out frequency scale
        freqRange = MRS_struct.p.sw(ii)/MRS_struct.p.LarmorFreq(ii);
        if MRS_struct.p.phantom
            F0 = 4.8;
        else
            F0 = 4.68;
        end
        MRS_struct.spec.freq = (MRS_struct.p.ZeroFillTo(ii) + 1 - (1:1:MRS_struct.p.ZeroFillTo(ii))) / MRS_struct.p.ZeroFillTo(ii) * freqRange + F0 - freqRange/2;
        
        MRS_struct.p.df(ii) = abs(MRS_struct.spec.freq(1) - MRS_struct.spec.freq(2));
        MRS_struct.p.SpecRes(ii) = MRS_struct.p.sw(ii)/MRS_struct.p.npoints(ii);
        MRS_struct.p.SpecResNominal(ii) = MRS_struct.p.sw(ii)/MRS_struct.p.ZeroFillTo(ii);
        MRS_struct.p.Tacq(ii) = 1/MRS_struct.p.SpecRes(ii);
        
        % Frame-by-frame determination of frequency of residual water and Cr (if HERMES or GSH editing)
        if MRS_struct.p.HERMES || any(strcmp(MRS_struct.p.target,'GSH'))
            F0freqRange = MRS_struct.spec.freq - 3.02 >= -0.1 & MRS_struct.spec.freq - 3.02 <= 0.1;
        else
            F0freqRange = MRS_struct.spec.freq - F0 >= -0.2 & MRS_struct.spec.freq - F0 <= 0.2;
        end
        [~,FrameMaxPos] = max(abs(real(AllFramesFT(F0freqRange,:))),[],1);
        F0freqRange = MRS_struct.spec.freq(F0freqRange);
        MRS_struct.spec.F0freq(ii,:) = F0freqRange(FrameMaxPos);
        
        % Estimate average amount of F0 offset
        if MRS_struct.p.HERMES || any(strcmp(MRS_struct.p.target,'GSH'))
            MRS_struct.out.AvgDeltaF0(ii) = mean(F0freqRange(FrameMaxPos) - 3.02);
        elseif any(strcmp(MRS_struct.p.vendor,{'Siemens_rda','Siemens_twix','Siemens_dicom'}))
            MRS_struct.out.AvgDeltaF0(ii) = mean(F0freqRange(FrameMaxPos) - 4.7); % Siemens assumes 4.7 ppm as F0
        else
            MRS_struct.out.AvgDeltaF0(ii) = mean(F0freqRange(FrameMaxPos) - F0);
        end
        
        % Frame-by-frame alignment
        switch MRS_struct.p.AlignTo
            case {'Cr','Cho','NAA'}
                [AllFramesFTrealign, MRS_struct] = AlignUsingPeak(AllFramesFT,MRS_struct);
            case 'H2O'
                [AllFramesFTrealign, MRS_struct] = AlignUsingH2O(AllFramesFT,MRS_struct);
            case 'SpecReg'
                [AllFramesFTrealign, MRS_struct] = Spectral_Registration(MRS_struct,0);
            case 'SpecRegDual'
                %Dual-channel Spectral Registration is applied separately to ON and OFF and they are coregistered after...
                [AllFramesFTrealign, MRS_struct] = Spectral_Registration(MRS_struct,0,1);
            case 'SpecRegHERMES'
                [AllFramesFTrealign, MRS_struct] = Spectral_Registration_HERMES(MRS_struct);
            case 'RobustSpecReg'
                [AllFramesFTrealign, MRS_struct] = Robust_Spectral_Registration(MRS_struct);
            case 'none' % GO (180224)
                % do nothing
                AllFramesFTrealign = AllFramesFT;
                MRS_struct.out.reject(:,ii) = zeros(1,size(AllFramesFT,2));
            otherwise
                error('AlignTo parameter in GannetPreInitialise.m not recognized. Check spelling.');
        end
        
        MRS_struct.spec.AllFramesFTrealign = AllFramesFTrealign;
        
        % Separate ON/OFF data and generate DIFF spectra
        
        if strcmp(MRS_struct.p.AlignTo,'RobustSpecReg')
            fprintf('\nPerforming weighted averaging of subspectra...\n');
        end
        
        if MRS_struct.p.HERMES
            
            MRS_struct.spec.(vox{kk}).subspec.A(ii,:) = mean(AllFramesFTrealign(:,1:4:size(AllFramesFTrealign,2)),2);
            MRS_struct.spec.(vox{kk}).subspec.B(ii,:) = mean(AllFramesFTrealign(:,2:4:size(AllFramesFTrealign,2)),2);
            MRS_struct.spec.(vox{kk}).subspec.C(ii,:) = mean(AllFramesFTrealign(:,3:4:size(AllFramesFTrealign,2)),2);
            MRS_struct.spec.(vox{kk}).subspec.D(ii,:) = mean(AllFramesFTrealign(:,4:4:size(AllFramesFTrealign,2)),2);
            
            for jj = 1:length(MRS_struct.p.target)
                
                if strcmp(MRS_struct.p.AlignTo,'RobustSpecReg') % determine weights for weighted averaging (MM: 190423)
                    
                    ON_inds  = find(MRS_struct.fids.ON_OFF(jj,:) == 1);
                    OFF_inds = find(MRS_struct.fids.ON_OFF(jj,:) == 0);
                    DIFFs = zeros(size(AllFramesFTrealign,1), size(AllFramesFTrealign,2)/4);
                    inds = 1:2;
                    
                    for ll = 1:size(DIFFs,2)
                        tmpON  = sum(AllFramesFTrealign(:,ON_inds(inds)),2);
                        tmpOFF = sum(AllFramesFTrealign(:,OFF_inds(inds)),2);
                        DIFFs(:,ll) = tmpON - tmpOFF;
                        inds = inds + 2;
                    end
                    
                    DIFFs = ifft(ifftshift(DIFFs,1),[],1);
                    DIFFs = fftshift(fft(DIFFs(1:MRS_struct.p.npoints(ii),:),[],1),1);
                    
                    freq = (size(DIFFs,1) + 1 - (1:size(DIFFs,1))) / size(DIFFs,1) * freqRange + 4.68 - freqRange/2;
                    freqLim = freq <= 4.25 & freq >= 1.8;
                    D = zeros(size(DIFFs,2));
                    
                    for ll = 1:size(DIFFs,2)
                        for mm = 1:size(DIFFs,2)
                            tmp = sum((real(DIFFs(freqLim,ll)) - real(DIFFs(freqLim,mm))).^2) / sum(freqLim);
                            if tmp == 0
                                D(ll,mm) = NaN;
                            else
                                D(ll,mm) = tmp;
                            end
                        end
                    end
                    
                    d = nanmean(D);
                    w = 1./d.^2;
                    w = repelem(w,2);
                    w = w/sum(w);
                    w = repmat(w, [size(AllFramesFTrealign,1) 1]);
                    
                    OFF = sum(w .* AllFramesFTrealign(:,MRS_struct.fids.ON_OFF(jj,:)==0),2);
                    ON  = sum(w .* AllFramesFTrealign(:,MRS_struct.fids.ON_OFF(jj,:)==1),2);
                    
                    MRS_struct.out.reject(:,ii) = zeros(size(AllFramesFTrealign,2),1);
                    
                else
                    
                    OFF = mean(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF(jj,:)==0)' & MRS_struct.out.reject(:,ii)==0),2);
                    ON  = mean(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF(jj,:)==1)' & MRS_struct.out.reject(:,ii)==0),2);
                    
                end
                
                MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).off(ii,:) = OFF;
                MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).on(ii,:)  = ON;
                
                MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:) = (ON - OFF)/2;
                MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:) = (mean(AllFramesFT(:,MRS_struct.fids.ON_OFF(jj,:)==1),2) - mean(AllFramesFT(:,MRS_struct.fids.ON_OFF(jj,:)==0),2))/2;
                
                % Edit-OFF,-OFF spectrum (for Cr referencing)
                if strcmp(MRS_struct.p.AlignTo,'RobustSpecReg') % determine weights for weighted averaging (MM: 190423)
                    
                    OFF_OFF_inds = find(all(MRS_struct.fids.ON_OFF'==0,2));
                    OFF_OFFs = ifft(ifftshift(AllFramesFTrealign(:,OFF_OFF_inds),1),[],1);
                    OFF_OFFs = fftshift(fft(OFF_OFFs(1:MRS_struct.p.npoints(ii),:),[],1),1);
                    D = zeros(size(OFF_OFFs,2));
                    
                    for ll = 1:size(OFF_OFFs,2)
                        for mm = 1:size(OFF_OFFs,2)
                            tmp = sum((real(OFF_OFFs(freqLim,ll)) - real(OFF_OFFs(freqLim,mm))).^2) / sum(freqLim);
                            if tmp == 0
                                D(ll,mm) = NaN;
                            else
                                D(ll,mm) = tmp;
                            end
                        end
                    end
                    
                    d = nanmean(D);
                    w = 1./d.^2;
                    w = w/sum(w);
                    w = repmat(w, [size(AllFramesFTrealign,1) 1]);
                    OFF_OFF = sum(w .* AllFramesFTrealign(:,OFF_OFF_inds),2);
                    
                else
                    
                    OFF_OFF = mean(AllFramesFTrealign(:,all(MRS_struct.fids.ON_OFF'==0,2) & MRS_struct.out.reject(:,ii)==0),2);
                    
                end
                
                MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).off_off(ii,:) = OFF_OFF;
                
            end
            
        else
            
            if strcmp(MRS_struct.p.AlignTo,'RobustSpecReg')
                
                % Determine weights for weighted averaging (MM: 190423)
                ON_inds  = find(MRS_struct.fids.ON_OFF == 1);
                OFF_inds = find(MRS_struct.fids.ON_OFF == 0);
                DIFFs = zeros(size(AllFramesFTrealign,1), size(AllFramesFTrealign,2)/2);
                
                for ll = 1:size(AllFramesFTrealign,2)/2
                    DIFFs(:,ll) = AllFramesFTrealign(:,ON_inds(ll)) - AllFramesFTrealign(:,OFF_inds(ll));
                end
                
                DIFFs = ifft(ifftshift(DIFFs,1),[],1);
                DIFFs = fftshift(fft(DIFFs(1:MRS_struct.p.npoints(ii),:),[],1),1);
                
                freq = (size(DIFFs,1) + 1 - (1:size(DIFFs,1))) / size(DIFFs,1) * freqRange + 4.68 - freqRange/2;
                freqLim = freq <= 4.25 & freq >= 1.8;
                D = zeros(size(AllFramesFTrealign,2)/2);
                
                for ll = 1:size(AllFramesFTrealign,2)/2
                    for mm = 1:size(AllFramesFTrealign,2)/2
                        tmp = sum((real(DIFFs(freqLim,ll)) - real(DIFFs(freqLim,mm))).^2) / sum(freqLim);
                        if tmp == 0
                            D(ll,mm) = NaN;
                        else
                            D(ll,mm) = tmp;
                        end
                    end
                end
                
                d = nanmean(D);
                w = 1./d.^2;
                w = w/sum(w);
                w = repmat(w, [size(AllFramesFTrealign,1) 1]);
                
                MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).off(ii,:) = sum(w .* AllFramesFTrealign(:,MRS_struct.fids.ON_OFF==0),2);
                MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).on(ii,:)  = sum(w .* AllFramesFTrealign(:,MRS_struct.fids.ON_OFF==1),2);
                
                MRS_struct.out.reject(:,ii) = zeros(size(AllFramesFTrealign,2),1);
                
            else
                
                MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).off(ii,:) = mean(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==0)' & MRS_struct.out.reject(:,ii)==0), 2);
                MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).on(ii,:)  = mean(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==1)' & MRS_struct.out.reject(:,ii)==0), 2);
                
            end
            
            MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).diff(ii,:) = (MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).on(ii,:) - MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).off(ii,:))/2;
            MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).diff_noalign(ii,:) = (mean(AllFramesFT(:,(MRS_struct.fids.ON_OFF==1)),2) - mean(AllFramesFT(:,(MRS_struct.fids.ON_OFF==0)),2))/2;
            
        end
        
        % Remove residual water from diff and diff_noalign spectra using HSVD -- GO & MGSaleh 2016
        if MRS_struct.p.water_removal
            
            for jj = 1:length(MRS_struct.p.target)
                if jj == 1
                    fprintf('\nFiltering out residual water signal...\n');
                end
                
                % Convert DIFF spectra to time domain, apply water filter, convert back to frequency domain
                MRS_struct.fids.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:) = waterremovalSVD(ifft(ifftshift(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:).')), ...
                    MRS_struct.p.sw(ii)/1e3, 8, -0.08, 0.08, 0, 2048);
                MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:) = fftshift(fft(MRS_struct.fids.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:)));
                
                MRS_struct.fids.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:) = waterremovalSVD(ifft(ifftshift(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:).')), ...
                    MRS_struct.p.sw(ii)/1e3, 8, -0.08, 0.08, 0, 2048);
                MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:) = fftshift(fft(MRS_struct.fids.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:)));
                
                % MM (170703): Need to perform baseline correction on filtered data
                freqbounds = MRS_struct.spec.freq <= 10 & MRS_struct.spec.freq >= 9;
                baseMean_diff = mean(real(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,freqbounds)));
                baseMean_diffnoalign = mean(real(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,freqbounds)));
                
                MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:) = MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:) - baseMean_diff;
                MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:) = MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:) - baseMean_diffnoalign;
            end
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   6. Build GannetLoad Output
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if ishandle(101)
            clf(101);
        end
        h = figure(101);
        % Open figure in center of screen
        scr_sz = get(0, 'ScreenSize');
        fig_w = 1000;
        fig_h = 707;
        set(h,'Position',[(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
        set(h,'Color',[1 1 1]);
        figTitle = 'GannetLoad Output';
        set(h,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');
        
        % Top left
        if length(MRS_struct.p.target) == 3
            subplot(5,2,1:2:5);
        else
            subplot(2,2,1);
        end
        GannetPlotPrePostAlign(MRS_struct, vox, ii, kk);               
        
        % Top right
        if MRS_struct.p.phantom
            F0 = 4.8;
        elseif MRS_struct.p.HERMES || any(strcmp(MRS_struct.p.target,'GSH'))
            F0 = 3.02;
        else
            F0 = 4.68;
        end
        
        subplot(2,2,2);
        rejectframesplot = (1./MRS_struct.out.reject(:,ii).') .* MRS_struct.spec.F0freq(ii,:);
        hold on;
        plot([1 size(FullData,2)], [F0 F0], '-k')
        plot([1 size(FullData,2)], [F0-0.04 F0-0.04], '--k')
        plot([1 size(FullData,2)], [F0+0.04 F0+0.04], '--k');
        plot(1:size(FullData,2), MRS_struct.spec.F0freq(ii,:)', '-', 1:size(FullData,2), rejectframesplot, 'ro');
        hold off;
        if MRS_struct.p.HERMES || any(strcmp(MRS_struct.p.target,'GSH'))
            text(size(FullData,2) + 0.025*size(FullData,2), F0, {'Nominal','Cr freq.'}, 'FontSize', 8);
        else
            text(size(FullData,2) + 0.025*size(FullData,2), F0, {'Nominal','water freq.'}, 'FontSize', 8);
        end
        set(gca,'TickDir','out','box','off','XLim',[1 size(FullData,2)], ...
            'YLim',[min([F0-0.06 MRS_struct.spec.F0freq(ii,:)-0.005]) max([F0+0.06 MRS_struct.spec.F0freq(ii,:)+0.005])]);
        if size(FullData,2) == 2
            set(gca,'XTick',[1 2]);
        end
        xlabel('average');
        ylabel('ppm');
        if MRS_struct.p.HERMES || any(strcmp(MRS_struct.p.target,'GSH'))
            title('Cr Frequency');
        else
            title('Water Frequency');
        end
        
        % Bottom left
        if length(MRS_struct.p.target) == 3
            subplot(3,2,5);
        else
            subplot(2,2,3);
        end
        if ~strcmp(MRS_struct.p.AlignTo,'no')
            CrFitLimLow = 2.72;
            CrFitLimHigh = 3.12;
            plotrange = MRS_struct.spec.freq <= CrFitLimHigh & MRS_struct.spec.freq >= CrFitLimLow;
            CrFitRange = sum(plotrange);
            plotrealign = [real(AllFramesFT(plotrange,:)); real(AllFramesFTrealign(plotrange,:))];
            % Don't display rejects
            plotrealign(CrFitRange+1:end,(MRS_struct.out.reject(:,ii).'==1)) = min(plotrealign(:));
            imagesc(plotrealign);
            title({'Cr Frequency','(pre- and post-alignment)'});
            xlabel('average');
            set(gca,'YTick',[1 CrFitRange CrFitRange+CrFitRange*(CrFitLimHigh-3.02)/(CrFitLimHigh-CrFitLimLow) CrFitRange*2], ...
                'YTickLabel', [CrFitLimHigh CrFitLimLow 3.02 CrFitLimLow], 'XLim', [1 size(FullData,2)], 'YLim', [1 CrFitRange*2], ...
                'TickDir','out','box','off');
            if size(FullData,2) == 2
                set(gca,'XTick',[1 2]);
            end
            % Add in labels for pre/post
            text(size(plotrealign,2)/18*17, 0.4*size(plotrealign,1), 'PRE', 'Color', [1 1 1], 'HorizontalAlignment', 'right');
            text(size(plotrealign,2)/18*17, 0.9*size(plotrealign,1), 'POST', 'Color', [1 1 1], 'HorizontalAlignment', 'right');
        else
            tmp = 'No realignment';
            text(0, 0.9, tmp, 'FontName', 'Arial');
        end
        
        % Bottom right
        subplot(2,2,4);
        axis off;
                
        % MM (180112)
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
            [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii*2-1});
        else
            [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii});
        end
        fname = [tmp tmp2];
        if length(fname) > 30
            fname = [fname(1:30) '...'];
        end
        
        text(0.25, 0.9, 'Filename: ', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
        text(0.275, 0.9, fname, 'FontName', 'Arial', 'FontSize', 13, 'Interpreter', 'none');
        
        tmp = [num2str(MRS_struct.p.Navg(ii)) ' averages'];
        text(0.25, 0.8, 'N_{avg}: ', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
        text(0.275, 0.8, tmp, 'FontName', 'Arial', 'FontSize', 13);
        
        tmp = [num2str(prod(MRS_struct.p.voxdim(ii,:))/1e3) ' mL'];
        text(0.25, 0.7, 'Volume: ', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
        text(0.275, 0.7, tmp, 'FontName', 'Arial', 'FontSize', 13);
        
        text(0.25, 0.6, 'Alignment: ', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
        text(0.275, 0.6, MRS_struct.p.AlignTo, 'FontName', 'Arial', 'FontSize', 13);
        
        tmp = [num2str(MRS_struct.p.LB) ' Hz'];
        text(0.25, 0.5, 'LB: ', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
        text(0.275, 0.5, tmp, 'FontName', 'Arial', 'FontSize', 13);
        
        text(0.25, 0.4, 'Rejects: ', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
        if strcmp(MRS_struct.p.AlignTo,'RobustSpecReg')
            text(0.275, 0.4, 'n/a - Weighted averaging used', 'FontName', 'Arial', 'FontSize', 13);
        else
            text(0.275, 0.4, num2str(sum(MRS_struct.out.reject(:,ii),1)), 'FontName', 'Arial', 'FontSize', 13);
        end
        
        text(0.25, 0.3, 'LoadVer: ', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
        text(0.275, 0.3, MRS_struct.version.load, 'FontName', 'Arial', 'FontSize', 13);
        
        % Gannet logo
        Gannet_path = which('GannetLoad');
        Gannet_logo = [Gannet_path(1:end-13) '/Gannet3_logo.png'];
        I = imread(Gannet_logo,'png','BackgroundColor',[1 1 1]);
        axes('Position',[0.80, 0.05, 0.15, 0.15]);
        imshow(I);
        text(0.9, 0, MRS_struct.version.Gannet, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
        axis off;
        axis square;
        
        % For Philips .data
        if strcmpi(MRS_struct.p.vendor,'Philips_data')
            fullpath = MRS_struct.metabfile{ii};
            fullpath = regexprep(fullpath, '.data', '_data');
            fullpath = regexprep(fullpath, '\', '_');
            fullpath = regexprep(fullpath, '/', '_');
        end
        
        % MM (180112)
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
            [~,metabfile_nopath] = fileparts(MRS_struct.metabfile{ii*2-1});
        else
            [~,metabfile_nopath] = fileparts(MRS_struct.metabfile{ii});
        end
        
        if any(strcmp(listfonts,'Arial'))
            set(findall(h,'-property','FontName'),'FontName','Arial');
        end
        
        % Create output folder
        if ~exist(fullfile(pwd, 'GannetLoad_output'),'dir')
            mkdir(fullfile(pwd, 'GannetLoad_output'));
        end
        
        % Save PDF
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[11 8.5]);
        set(h,'PaperPosition',[0 0 11 8.5]);
        
        if strcmpi(MRS_struct.p.vendor,'Philips_data')
            pdfname = fullfile(pwd, 'GannetLoad_output', [fullpath '_' vox{kk} '_load.pdf']);
        else
            pdfname = fullfile(pwd, 'GannetLoad_output', [metabfile_nopath '_' vox{kk} '_load.pdf']);
        end
        saveas(h, pdfname);
        
        % Export the processed data into an SDAT file
        if MRS_struct.p.sdat
            if strcmpi(MRS_struct.p.vendor,'Philips')
                % Set up filenames
                sdat_G_name = ['GannetLoad_output/' metabfile_nopath  '_' vox{kk} '_G.sdat'];
                spar_G_name = ['GannetLoad_output/' metabfile_nopath  '_' vox{kk} '_G.spar'];
                % Make file copies for SDAT/SPAR files
                copyfile(metabfile{ii},sdat_G_name);
                sparname = [metabfile{ii}(1:end-4) MRS_struct.p.spar_string];
                copyfile(sparname,spar_G_name);
                % Write DIFF data into the SDAT file
                sdat_diff_out = conj(ifft(fftshift(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).diff(ii,:),2),[],2));
                sdat_diff_out = sdat_diff_out(1:MRS_struct.p.npoints(ii));
                % Also write out OFF data
                sdat_off_out = conj(ifft(fftshift(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).off(ii,:),2),[],2));
                sdat_off_out = sdat_off_out(1:MRS_struct.p.npoints(ii));
                fileid  = fopen(sdat_G_name,'w','ieee-le');
                ff(:,1:2:2*MRS_struct.p.npoints(ii)) = real(sdat_diff_out);
                ff(:,2:2:2*MRS_struct.p.npoints(ii)) = imag(sdat_diff_out);
                gg(:,1:2:2*MRS_struct.p.npoints(ii)) = real(sdat_off_out);
                gg(:,2:2:2*MRS_struct.p.npoints(ii)) = imag(sdat_off_out);
                fwriteVAXD(fileid,[ff.' gg.'],'float');
                fclose(fileid);
            else
                warning('Only Philips SDAT files can be exported! No data exported.');
            end      
        end
        
        % 140116: ADH reorder structure
        if isfield(MRS_struct, 'waterfile')
            structorder = {'version', 'ii', 'metabfile', ...
                'waterfile', 'p', 'fids', 'spec', 'out'};
        else
            structorder = {'version', 'ii', 'metabfile', ...
                'p', 'fids', 'spec', 'out'};
        end
        MRS_struct = orderfields(MRS_struct, structorder);
        
        % Save MRS_struct as mat file
        if ii == numscans && MRS_struct.p.mat
            % Set up filename
            mat_name = ['GannetLoad_output/MRS_struct_' vox{kk} '.mat'];
            save(mat_name,'MRS_struct');
        end
        
    end % end of output loop over voxels
    
end % end of load-and-processing loop over datasets

end



