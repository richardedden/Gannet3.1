function MRS_struct = GannetFit(MRS_struct, varargin)
% Gannet 3.1 GannetFit
% Started by RAEE Nov 5, 2012
% Updates by MGS, MM 2016-2019

MRS_struct.version.fit = '190530';

if MRS_struct.p.PRIAM
    vox = MRS_struct.p.Vox;
else
    vox = MRS_struct.p.Vox(1);
end

if nargin < 2
    target = MRS_struct.p.target;
elseif nargin > 1
    % varargin = Optional arguments if user wants to overwrite fitting
    %            parameters set in GannetPreInitialise; can include several
    %            options, which are:
    %            'GABA' or 'Glx': target metabolite
    switch varargin{1}
        case 'GABA'
            MRS_struct.p.target = 'GABA';
        case 'Glx'
            MRS_struct.p.target = 'Glx';
        case 'GABAGlx'
            MRS_struct.p.target = 'GABAGlx';
        case 'GSH'
            MRS_struct.p.target = 'GSH';
        case 'Lac'
            MRS_struct.p.target = 'Lac';
        case 'EtOH'
            MRS_struct.p.target = 'EtOH';
    end
    target = {MRS_struct.p.target};
end

freq = MRS_struct.spec.freq;

lsqopts = optimset('lsqcurvefit');
lsqopts = optimset(lsqopts,'MaxIter',800,'TolX',1e-4,'TolFun',1e-4,'Display','off');
nlinopts = statset('nlinfit');
nlinopts = statset(nlinopts,'MaxIter',400,'TolX',1e-6,'TolFun',1e-6,'FunValCheck','off');

warning('on');
warning('off','stats:nlinfit:ModelConstantWRTParam');
warning('off','stats:nlinfit:IllConditionedJacobian');
warning('off','MATLAB:rankDeficientMatrix');

% Loop over voxels if PRIAM
for kk = 1:length(vox)
    
    if strcmp(MRS_struct.p.Reference_compound,'H2O')
        WaterData = MRS_struct.spec.(vox{kk}).water;
    end
    
    % Loop over edited spectra if HERMES
    for jj = 1:length(target)
        
        fprintf('\nFitting %s...',target{jj});
        
        DIFF = MRS_struct.spec.(vox{kk}).(target{jj}).diff;
        if MRS_struct.p.HERMES
            OFF = MRS_struct.spec.(vox{kk}).(target{jj}).off_off;
        else
            OFF = MRS_struct.spec.(vox{kk}).(target{jj}).off;
        end
        
        numscans = MRS_struct.p.numscans;
        
        for ii = 1:numscans
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   1.  Metabolite Fitting
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if strcmp(target{jj},'GABA')
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   GABA
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if length(MRS_struct.p.target) == 3 && all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
                    
                    freqbounds = find(freq <= 3.55 & freq >= 2.6);
                    plotbounds = find(freq <= 3.6 & freq >= 2.6);
                    
                    maxinGABA = abs(max(real(DIFF(ii,freqbounds))) - min(real(DIFF(ii,freqbounds))));
                    grad_points = (real(DIFF(ii,freqbounds(end))) - real(DIFF(ii,freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1));
                    LinearInit = grad_points ./ abs(freq(1) - freq(2));
                    constInit = (real(DIFF(ii,freqbounds(end))) + real(DIFF(ii,freqbounds(1))))./2;
                    
                    GaussModelInit = [maxinGABA -90 3.026 -LinearInit constInit];
                    GaussModelInit([1 4 5]) = GaussModelInit([1 4 5]) / maxinGABA; % Scale initial conditions to avoid warnings about numerical underflow
                    
                    lb = [0 -200 2.87 -40*maxinGABA -2000*maxinGABA];
                    ub = [4000*maxinGABA -40 3.12 40*maxinGABA 1000*maxinGABA];
                    lb([1 4 5]) = lb([1 4 5]) / maxinGABA;
                    ub([1 4 5]) = ub([1 4 5]) / maxinGABA;
                    
                    % Down-weight co-edited Cho signal by including
                    % observation weights in nonlinear regression
                    w = ones(size(DIFF(ii,freqbounds)));
                    residfreq = freq(freqbounds);
                    ChoRange = residfreq >= 3.16 & residfreq <= 3.285;
                    weightRange = ChoRange;
                    w(weightRange) = 0.001;
                    
                    % Least-squares model fitting
                    GaussModelInit = lsqcurvefit(@GaussModel, GaussModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinGABA, lb, ub, lsqopts);
                    modelFun_w = @(x,freq) sqrt(w) .* GaussModel(x,freq); % add weights to the model
                    [GaussModelParam, resid] = nlinfit(freq(freqbounds), sqrt(w) .* real(DIFF(ii,freqbounds)) / maxinGABA, modelFun_w, GaussModelInit, nlinopts); % add weights to the data
                    [~, residPlot] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinGABA, @GaussModel, GaussModelParam, nlinopts); % re-run for residuals for output figure
                    
                    % Rescale fit parameters and residuals
                    GaussModelParam([1 4 5]) = GaussModelParam([1 4 5]) * maxinGABA;
                    resid = resid * maxinGABA;
                    residPlot = residPlot * maxinGABA;
                    
                else
                    
                    freqbounds = find(freq <= 3.55 & freq >= 2.79);
                    plotbounds = find(freq <= 3.6 & freq >= 2.7);
                    
                    maxinGABA = abs(max(real(DIFF(ii,freqbounds))) - min(real(DIFF(ii,freqbounds))));
                    grad_points = (real(DIFF(ii,freqbounds(end))) - real(DIFF(ii,freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1));
                    LinearInit = grad_points ./ abs(freq(1) - freq(2));
                    constInit = (real(DIFF(ii,freqbounds(end))) + real(DIFF(ii,freqbounds(1))))./2;
                    
                    GaussModelInit = [maxinGABA -90 3.026 -LinearInit constInit];
                    GaussModelInit([1 4 5]) = GaussModelInit([1 4 5]) / maxinGABA; % Scale initial conditions to avoid warnings about numerical underflow
                    
                    lb = [0 -200 2.87 -40*maxinGABA -2000*maxinGABA];
                    ub = [4000*maxinGABA -40 3.12 40*maxinGABA 1000*maxinGABA];
                    lb([1 4 5]) = lb([1 4 5]) / maxinGABA;
                    ub([1 4 5]) = ub([1 4 5]) / maxinGABA;
                    
                    % Down-weight Cho subtraction artifact by including
                    % observation weights in nonlinear regression; improves
                    % accuracy of peak fitting (MM: 170701 - thanks to Alex
                    % Craven of University of Bergen for this idea)
                    w = ones(size(DIFF(ii,freqbounds)));
                    residfreq = freq(freqbounds);
                    ChoRange = residfreq >= 3.16 & residfreq <= 3.285;
                    weightRange = ChoRange;
                    w(weightRange) = 0.001;
                    
                    % Weighted least-squares model fitting
                    GaussModelInit = lsqcurvefit(@GaussModel, GaussModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinGABA, lb, ub, lsqopts);
                    modelFun_w = @(x,freq) sqrt(w) .* GaussModel(x,freq); % add weights to the model
                    [GaussModelParam, resid] = nlinfit(freq(freqbounds), sqrt(w) .* real(DIFF(ii,freqbounds)) / maxinGABA, modelFun_w, GaussModelInit, nlinopts); % add weights to the data
                    [~, residPlot] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinGABA, @GaussModel, GaussModelParam, nlinopts); % re-run for residuals for output figure
                    
                    % Rescale fit parameters and residuals
                    GaussModelParam([1 4 5]) = GaussModelParam([1 4 5]) * maxinGABA;
                    resid = resid * maxinGABA;
                    residPlot = residPlot * maxinGABA;
                    
                end
                
                GABAheight = GaussModelParam(1);
                MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii) = 100*std(resid)/GABAheight;
                MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) = GaussModelParam(1)./sqrt(-GaussModelParam(2))*sqrt(pi);
                sigma = sqrt(1/(2*(abs(GaussModelParam(2)))));
                MRS_struct.out.(vox{kk}).(target{jj}).FWHM(ii) = abs((2*MRS_struct.p.LarmorFreq(ii))*sigma);
                MRS_struct.out.(vox{kk}).(target{jj}).ModelParam(ii,:) = GaussModelParam;
                MRS_struct.out.(vox{kk}).(target{jj}).Resid(ii,:) = resid;
                
                % Calculate SNR of GABA signal
                noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
                MRS_struct.out.(vox{kk}).(target{jj}).SNR(ii) = abs(GABAheight)/noiseSigma_DIFF;
                
            elseif strcmp(target{jj},'GSH')
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   GSH
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                freqbounds = find(freq <= 3.5 & freq >= 2.25);
                plotbounds = find(freq <= 4.2 & freq >= 1.75);
                
                GSHbounds = freq <= 3.3 & freq >= 2.85;
                Aspartylbounds = freq <= 2.85 & freq >= 2.25;
                
                maxinGSH = max(abs(real(DIFF(ii,GSHbounds))));
                [maxinAspartyl, maxInd] = max(abs(real(DIFF(ii,Aspartylbounds))));
                
                offset = real(DIFF(ii,freqbounds(end)));
                grad_points = (real(DIFF(ii,freqbounds(end))) - real(DIFF(ii,freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1));
                LinearInit = grad_points ./ abs(freq(1) - freq(2));
                
                tmp = DIFF(ii,Aspartylbounds);
                s = sign(real(tmp(maxInd)));
                maxinAspartyl = s * maxinAspartyl;
                
                if MRS_struct.p.HERMES
                    s = -1;
                else
                    s = 1;
                end
                
                if strcmp(MRS_struct.p.GSH_model,'FiveGauss')
                    
                    GSHgaussModel = @FiveGaussModel;
                    
                    GaussModelInit = [maxinGSH             -300  2.95 ...
                                      s*maxinAspartyl*0.25 -500  2.73 ...
                                      maxinAspartyl        -1000 2.61 ...
                                      maxinAspartyl        -1000 2.55 ...
                                      s*maxinAspartyl*0.15 -600  2.45 ...
                                      offset -LinearInit -LinearInit];
                    GaussModelInit([1 4 7 10 13 16 17 18]) = GaussModelInit([1 4 7 10 13 16 17 18]) / maxinGSH; % Scale initial conditions to avoid warnings about numerical underflow
                    
                    lb = [-4000*maxinGSH            -1000 2.95-0.02 ...
                           4000*maxinAspartyl*0.25  -1000 2.73-0.02 ...
                           4000*maxinAspartyl       -1000 2.61-0.02 ...
                           4000*maxinAspartyl       -1000 2.55-0.02 ...
                           4000*maxinAspartyl*0.15  -1000 2.45-0.02 ...
                          -2000*abs(offset) 2000*maxinAspartyl 2000*maxinAspartyl];
                    ub =  [4000*maxinGSH           -40 2.95+0.02 ...
                          -4000*maxinAspartyl*0.25 -40 2.73+0.02 ...
                          -4000*maxinAspartyl      -40 2.61+0.02 ...
                          -4000*maxinAspartyl      -40 2.55+0.02 ...
                          -4000*maxinAspartyl*0.15 -40 2.45+0.02 ...
                           1000*abs(offset) -1000*maxinAspartyl -1000*maxinAspartyl];
                    lb([1 4 7 10 13 16 17 18]) = lb([1 4 7 10 13 16 17 18]) / maxinGSH;
                    ub([1 4 7 10 13 16 17 18]) = ub([1 4 7 10 13 16 17 18]) / maxinGSH;
                    
                elseif strcmp(MRS_struct.p.GSH_model,'SixGauss')
                    
                    GSHgaussModel = @SixGaussModel;
                    
                    GaussModelInit = [maxinGSH           -300  2.95 ...
                                      maxinAspartyl*0.7  -500  2.73 ...
                                      maxinAspartyl      -1000 2.63 ...
                                      maxinAspartyl*0.7  -1000 2.58 ...
                                      maxinAspartyl*0.5  -600  2.46 ...
                                      maxinAspartyl*0.35 -600  2.37 ...
                                      offset -LinearInit -LinearInit];
                    GaussModelInit([1 4 7 10 13 16 19 20 21]) = GaussModelInit([1 4 7 10 13 16 19 20 21]) / maxinGSH; % Scale initial conditions to avoid warnings about numerical underflow
                    
                    lb = [-4000*maxinGSH           -1000 2.95-0.02 ...
                          -4000*maxinAspartyl*0.7  -1000 2.73-0.02 ...
                          -4000*maxinAspartyl      -1000 2.63-0.02 ...
                          -4000*maxinAspartyl*0.7  -1000 2.58-0.02 ...
                          -4000*maxinAspartyl*0.5  -1000 2.46-0.02 ...
                          -4000*maxinAspartyl*0.35 -1000 2.37-0.02 ...
                          -2000*abs(offset) -2000*maxinAspartyl -2000*maxinAspartyl];
                    ub =  [4000*maxinGSH           -40 2.95+0.02 ...
                           4000*maxinAspartyl*0.7  -40 2.73+0.02 ...
                           4000*maxinAspartyl      -40 2.63+0.02 ...
                           4000*maxinAspartyl*0.7  -40 2.58+0.02 ...
                           4000*maxinAspartyl*0.5  -40 2.46+0.02 ...
                           4000*maxinAspartyl*0.35 -40 2.37+0.02 ...
                           1000*abs(offset) 1000*maxinAspartyl 1000*maxinAspartyl];
                    lb([1 4 7 10 13 16 19 20 21]) = lb([1 4 7 10 13 16 19 20 21]) / maxinGSH;
                    ub([1 4 7 10 13 16 19 20 21]) = ub([1 4 7 10 13 16 19 20 21]) / maxinGSH;
                    
                end
                
                if length(MRS_struct.p.target) == 3 && all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
                    w = ones(size(DIFF(ii,freqbounds)));
                    residfreq = freq(freqbounds);
                    ChoRange = residfreq >= 3.13 & residfreq <= 3.3;
                    weightRange = ChoRange;
                    w(weightRange) = 0.001;
                else
                    w = ones(size(DIFF(ii,freqbounds)));
                end
                
                % Weighted least-squares model fitting
                GaussModelInit = lsqcurvefit(GSHgaussModel, GaussModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinGSH, lb, ub, lsqopts);
                modelFun_w = @(x,freq) sqrt(w) .* GSHgaussModel(x,freq); % add weights to the model
                [GaussModelParam, resid] = nlinfit(freq(freqbounds), sqrt(w) .* real(DIFF(ii,freqbounds)) / maxinGSH, modelFun_w, GaussModelInit, nlinopts); % add weights to the data
                [~, residPlot] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinGSH, GSHgaussModel, GaussModelParam, nlinopts); % re-run for residuals for output figure
                
                % Rescale fit parameters and residuals
                if strcmp(MRS_struct.p.GSH_model,'FiveGauss')
                    GaussModelParam([1 4 7 10 13 16 17 18]) = GaussModelParam([1 4 7 10 13 16 17 18]) * maxinGSH;
                elseif strcmp(MRS_struct.p.GSH_model,'SixGauss')
                    GaussModelParam([1 4 7 10 13 16 19 20 21]) = GaussModelParam([1 4 7 10 13 16 19 20 21]) * maxinGSH;
                end
                resid = resid * maxinGSH;
                residPlot = residPlot * maxinGSH;
                
                GSHGaussModelParam = GaussModelParam;
                if strcmp(MRS_struct.p.GSH_model,'FiveGauss')
                    GSHGaussModelParam(4:3:13) = 0;
                elseif strcmp(MRS_struct.p.GSH_model,'SixGauss')
                    GSHGaussModelParam(4:3:16) = 0;
                end
                
                BaselineModelParam = GSHGaussModelParam;
                BaselineModelParam(1) = 0;
                
                MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) = real(sum(FiveGaussModel(GSHGaussModelParam, freq(freqbounds)) - FiveGaussModel(BaselineModelParam, freq(freqbounds)))) ...
                    * abs(freq(1)-freq(2));
                GSHheight = GSHGaussModelParam(1);
                
                % Range to determine residuals for GSH
                residfreq = freq(freqbounds);
                residGSH = resid(residfreq <= 3.3 & residfreq >= 2.82);
                
                MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii) = 100*std(residGSH)/GSHheight;
                sigma = sqrt(1/(2*(abs(GSHGaussModelParam(2)))));
                MRS_struct.out.(vox{kk}).(target{jj}).FWHM(ii) =  abs((2*MRS_struct.p.LarmorFreq(ii))*sigma);
                MRS_struct.out.(vox{kk}).(target{jj}).ModelParam(ii,:) = GaussModelParam;
                MRS_struct.out.(vox{kk}).(target{jj}).Resid(ii,:) = residGSH;
                
                % Calculate SNR of GSH signal
                noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
                MRS_struct.out.(vox{kk}).(target{jj}).SNR(ii) = abs(GSHheight)/noiseSigma_DIFF;
                
            elseif strcmp(target{jj},'Lac')
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   Lac
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                freqbounds = find(freq <= 1.8 & freq >= 0.5);
                plotbounds = find(freq <= 2 & freq >= 0);
                
                offset = (mean(real(DIFF(ii,freqbounds(1:10))),2) + mean(real(DIFF(ii,freqbounds((end-9):end))),2))/2;
                slope = (mean(real(DIFF(ii,freqbounds(1:10))),2) - mean(real(DIFF(ii,freqbounds((end-9):end))),2))/abs(freq(freqbounds(1)) - freq(freqbounds(end)));
                peak_amp = 0.03; % Presumably this won't work for some data... for now it seems to work
                
                FourGaussModelInit = [peak_amp*0.16 -100 1.18 peak_amp*0.3 -1000 1.325 offset slope 0];
                lb = [0 -300 0.9 0 -5000 1.0  -1 -1 -1];
                ub = [1 0 1.4 1 0 1.6  1 1 1];
                
                FourGaussModelInit = lsqcurvefit(@FourGaussModel, FourGaussModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)), lb, ub,lsqopts);
                [FourGaussModelParam, resid] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)), @FourGaussModel, FourGaussModelInit, nlinopts);
                
                MRS_struct.out.(vox{kk}).Lac.ModelParam(ii,:) = FourGaussModelParam;
                %LacGaussModelParam = FourGaussModelParam;
                %LacGaussModelParam(1) = 0;
                MMGaussModelParam = FourGaussModelParam;
                MMGaussModelParam(4) = 0;
                
                MRS_struct.out.(vox{kk}).Lac.Area(ii) = real(sum(FourGaussModel([FourGaussModelParam(1:6) 0 0 0],freq(freqbounds)))) * abs(freq(1)-freq(2)); % NB: this is Lac+MM
                Lacheight = FourGaussModelParam(4);
                MRS_struct.out.(vox{kk}).Lac.FitError(ii) = 100*std(resid)/Lacheight;
                MRS_struct.out.(vox{kk}).Lac.FWHM(ii) = NaN; % MM (170818): Still need to calculate FWHM
                MRS_struct.out.(vox{kk}).Lac.Resid(ii,:) = resid;
                
                % Calculate SNR of Lac signal (MM: 170502)
                noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
                MRS_struct.out.(vox{kk}).Lac.SNR(ii) = abs(Lacheight)/noiseSigma_DIFF;
                
            elseif strcmp(target{jj},'Glx')
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   Glx
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                freqbounds = find(freq <= 4.1 & freq >= 3.45);
                plotbounds = find(freq <= 4.5 & freq >= 3);
                
                maxinGABA = max(real(DIFF(ii,freqbounds)));
                grad_points = (real(DIFF(ii,freqbounds(end))) - real(DIFF(ii,freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1));
                LinearInit = grad_points ./ abs(freq(1) - freq(2));
                constInit = (real(DIFF(ii,freqbounds(end))) + real(DIFF(ii,freqbounds(1))))./2;
                
                GaussModelInit = [maxinGABA -90 3.72 maxinGABA -90 3.77 -LinearInit constInit];
                lb = [0 -200 3.72-0.01 0 -200 3.77-0.01 -40*maxinGABA -2000*maxinGABA];
                ub = [4000*maxinGABA -40 3.72+0.01 4000*maxinGABA -40 3.77+0.01 40*maxinGABA 1000*maxinGABA];
                
                GaussModelInit = lsqcurvefit(@DoubleGaussModel, GaussModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)), lb, ub, lsqopts);
                [GaussModelParam, resid] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)), @DoubleGaussModel, GaussModelInit, nlinopts);
                
                Glxheight = max(GaussModelParam([1,4]));
                MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii) = 100*std(resid)/Glxheight;
                MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) = (GaussModelParam(1)./sqrt(-GaussModelParam(2))*sqrt(pi)) + ...
                    (GaussModelParam(4)./sqrt(-GaussModelParam(5))*sqrt(pi));
                sigma = ((1/(2*(abs(GaussModelParam(2))))).^(1/2)) + ((1/(2*(abs(GaussModelParam(5))))).^(1/2));
                MRS_struct.out.(vox{kk}).(target{jj}).FWHM(ii) = abs((2*MRS_struct.p.LarmorFreq(ii))*sigma);
                MRS_struct.out.(vox{kk}).(target{jj}).ModelParam(ii,:) = GaussModelParam;
                MRS_struct.out.(vox{kk}).(target{jj}).Resid(ii,:) = resid;
                
                % Calculate SNR of Glx signal
                noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
                MRS_struct.out.(vox{kk}).Glx.SNR(ii) = abs(Glxheight)/noiseSigma_DIFF;
                
            elseif strcmp(target{jj},'GABAGlx')
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   GABA+Glx
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                freqbounds = find(freq <= 4.1 & freq >= 2.79);
                plotbounds = find(freq <= 4.2 & freq >= 2.7);
                
                GABAbounds = freq <= 3.2 & freq >= 2.79;
                Glxbounds  = freq <= 4.1 & freq >= 3.4;
                
                maxinGABA = max(real(DIFF(ii,GABAbounds)));
                maxinGlx = max(real(DIFF(ii,Glxbounds)));
                grad_points = (real(DIFF(ii,freqbounds(end))) - real(DIFF(ii,freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1));
                LinearInit = grad_points ./ abs(freq(1) - freq(2));
                
                GaussModelInit = [maxinGlx -700 3.71 maxinGlx -700 3.79 maxinGABA -90 3.02 -LinearInit 0 0];
                GaussModelInit([1 4 7 10]) = GaussModelInit([1 4 7 10]) / maxinGlx; % Scale initial conditions to avoid warnings about numerical underflow
                
                lb = [-4000*maxinGlx -1000 3.71-0.02 -4000*maxinGlx -1000 3.79-0.02 -4000*maxinGABA -200 3.02-0.05 -40*maxinGABA -2000*maxinGABA -2000*maxinGABA];
                ub = [4000*maxinGlx -40 3.71+0.02 4000*maxinGlx -40 3.79+0.02 4000*maxinGABA -40 3.02+0.05 40*maxinGABA 1000*maxinGABA 1000*maxinGABA];
                lb([1 4 7 10]) = lb([1 4 7 10]) / maxinGlx;
                ub([1 4 7 10]) = ub([1 4 7 10]) / maxinGlx;
                
                % Down-weight Cho subtraction artifact and (if HERMES)
                % signals downfield of Glx by including observation weights
                % in nonlinear regression; improves accuracy of peak
                % fittings (MM: 170701 - thanks to Alex Craven of
                % University of Bergen for this idea)
                w = ones(size(DIFF(ii,freqbounds)));
                residfreq = freq(freqbounds);
                ChoRange = residfreq >= 3.16 & residfreq <= 3.285;
                GlxDownFieldRange = residfreq >= 3.9 & residfreq <= 4.2;
                if MRS_struct.p.HERMES && any(strcmp(MRS_struct.p.vendor,{'Philips','Philips_data','Philips_raw'}))
                    weightRange = ChoRange | GlxDownFieldRange;
                else
                    weightRange = ChoRange;
                end
                w(weightRange) = 0.001;
                
                % Weighted least-squares model fitting
                GaussModelInit = lsqcurvefit(@GABAGlxModel, GaussModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinGlx, lb, ub, lsqopts);
                modelFun_w = @(x,freq) sqrt(w) .* GABAGlxModel(x,freq); % add weights to the model
                [GaussModelParam, resid] = nlinfit(freq(freqbounds), sqrt(w) .* real(DIFF(ii,freqbounds)) / maxinGlx, modelFun_w, GaussModelInit, nlinopts); % add weights to the data
                [~, residPlot] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinGlx, @GABAGlxModel, GaussModelParam, nlinopts); % re-run for residuals for output figure
                
                % Rescale fit parameters and residuals
                GaussModelParam([1 4 7 10 11 12]) = GaussModelParam([1 4 7 10 11 12]) * maxinGlx;
                resid = resid * maxinGlx;
                residPlot = residPlot * maxinGlx;
                
                % Range to determine residuals for GABA and Glx
                residGABA = resid(residfreq <= 3.55 & residfreq >= 2.79);
                residGlx = resid(residfreq <= 4.10 & residfreq >= 3.45);
                
                % GABA fitting output
                MRS_struct.out.(vox{kk}).GABA.Area(ii) = (GaussModelParam(7)./sqrt(-GaussModelParam(8))*sqrt(pi));
                GABAheight = GaussModelParam(7);
                MRS_struct.out.(vox{kk}).GABA.FitError(ii) = 100*std(residGABA)/GABAheight;
                sigma = sqrt(1/(2*(abs(GaussModelParam(8)))));
                MRS_struct.out.(vox{kk}).GABA.FWHM(ii) = abs((2*MRS_struct.p.LarmorFreq(ii))*sigma);
                MRS_struct.out.(vox{kk}).GABA.ModelParam(ii,:) = GaussModelParam;
                MRS_struct.out.(vox{kk}).GABA.Resid(ii,:) = residGABA;
                
                % Calculate SNR of GABA signal
                noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
                MRS_struct.out.(vox{kk}).GABA.SNR(ii) = abs(GABAheight)/noiseSigma_DIFF;
                
                % Glx fitting output
                MRS_struct.out.(vox{kk}).Glx.Area(ii) = (GaussModelParam(1)./sqrt(-GaussModelParam(2))*sqrt(pi)) + ...
                    (GaussModelParam(4)./sqrt(-GaussModelParam(5))*sqrt(pi));
                Glxheight = max(GaussModelParam([1,4]));
                MRS_struct.out.(vox{kk}).Glx.FitError(ii) = 100*std(residGABA)/Glxheight;
                sigma = sqrt(1/(2*(abs(GaussModelParam(2))))) + sqrt(1/(2*(abs(GaussModelParam(5)))));
                MRS_struct.out.(vox{kk}).Glx.FWHM(ii) = abs((2*MRS_struct.p.LarmorFreq(ii))*sigma);
                MRS_struct.out.(vox{kk}).Glx.ModelParam(ii,:) = GaussModelParam;
                MRS_struct.out.(vox{kk}).Glx.Resid(ii,:) = residGlx;
                
                % Calculate SNR of Glx signal
                MRS_struct.out.(vox{kk}).Glx.SNR(ii) = abs(Glxheight)/noiseSigma_DIFF;
                
            elseif strcmp(target{jj},'EtOH')
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   EtOH
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                freqbounds = find(freq <= 1.8 & freq >= 0.6);
                plotbounds = find(freq <= 1.9 & freq >= 0.4);
                
                maxinEtOH = max(real(DIFF(ii,freqbounds)));
                grad_points = (real(DIFF(ii,freqbounds(end))) - real(DIFF(ii,freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1));
                LinearInit = grad_points ./ abs(freq(1) - freq(2));
                
                LorentzModelInit = [maxinEtOH 1.11 1/500 ...
                                    maxinEtOH 1.23 1/500 ...
                                    -LinearInit 0];
                LorentzModelInit([1 4 7]) = LorentzModelInit([1 4 7]) / maxinEtOH; % Scale initial conditions to avoid warnings about numerical underflow
                
                lb = [0             1.11-0.01 1/700 0             1.23-0.01 1/700 -40*maxinEtOH -2e3*maxinEtOH];
                ub = [100*maxinEtOH 1.11+0.01 1/300 100*maxinEtOH 1.23+0.01 1/300  40*maxinEtOH  1e3*maxinEtOH];
                
                lb([1 4 7]) = lb([1 4 7]) / maxinEtOH;
                ub([1 4 7]) = ub([1 4 7]) / maxinEtOH;
                
                % Down-weight co-edited Lac signal by including observation
                % weights in nonlinear regression
                w = ones(size(DIFF(ii,freqbounds)));
                residfreq = freq(freqbounds);
                LacRange = residfreq >= 1.29 & residfreq <= 1.51;
                weightRange = LacRange;
                w(weightRange) = 0.001;
                
                % Weighted least-squares model fitting
                LorentzModelInit = lsqcurvefit(@EtOHModel, LorentzModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinEtOH, lb, ub, lsqopts);
                modelFun_w = @(x,freq) sqrt(w) .* EtOHModel(x,freq); % add weights to the model
                [LorentzModelParam, resid] = nlinfit(freq(freqbounds), sqrt(w) .* real(DIFF(ii,freqbounds)) / maxinEtOH, modelFun_w, LorentzModelInit, nlinopts); % add weights to the data
                [~, residPlot] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinEtOH, @EtOHModel, LorentzModelParam, nlinopts); % re-run for residuals for output figure
                
                % Rescale fit parameters and residuals
                LorentzModelParam([1 4 7 8]) = LorentzModelParam([1 4 7 8]) * maxinEtOH;
                resid = resid * maxinEtOH;
                residPlot = residPlot * maxinEtOH;
                
                EtOHheight = max(EtOHModel([LorentzModelParam(1:end-2) 0 0],freq(freqbounds)));
                MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii) = 100*std(resid)/EtOHheight;
                MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) = sum(EtOHModel([LorentzModelParam(1:end-2) 0 0],freq(freqbounds))) * abs(freq(1)-freq(2));
                
                MRS_struct.out.(vox{kk}).(target{jj}).FWHM(ii) = (LorentzModelParam(3) + LorentzModelParam(6)) * MRS_struct.p.LarmorFreq(ii);
                MRS_struct.out.(vox{kk}).(target{jj}).ModelParam(ii,:) = LorentzModelParam;
                MRS_struct.out.(vox{kk}).(target{jj}).Resid(ii,:) = resid;
                
                noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
                MRS_struct.out.(vox{kk}).EtOH.SNR(ii) = abs(EtOHheight)/noiseSigma_DIFF;
                
            else
                
                error('Fitting %s not recognised',target{jj});
                
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   1a. Start up the output figure
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if ishandle(102)
                clf(102);
            end
            h = figure(102);
            % Open figure in center of screen
            scr_sz = get(0, 'ScreenSize');
            fig_w = 1000;
            fig_h = 707;
            set(h,'Position',[(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
            set(h,'Color',[1 1 1]);
            figTitle = 'GannetFit Output';
            set(h,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');
            
            % Spectra plot
            subplot(2,2,1);
            metabmin = min(real(DIFF(ii,plotbounds)));
            metabmax = max(real(DIFF(ii,plotbounds)));
            resmax = max(resid);
            resid = resid + metabmin - resmax;
            
            switch target{jj}
                case 'GABA'
                    residPlot = residPlot + metabmin - max(residPlot);
                    residPlot2 = residPlot;
                    residPlot2(weightRange) = NaN;
                    hold on;
                    plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b', ...
                        freq(freqbounds), GaussModel(GaussModelParam,freq(freqbounds)), 'r', ...
                        freq(freqbounds), residPlot2, 'k');
                    plot(freq(freqbounds(ChoRange)), residPlot(ChoRange), 'Color', [255 160 64]/255);
                    hold off;
                    set(gca,'XLim',[2.6 3.6]);
                    
                case 'GSH'
                    residPlot = residPlot + metabmin - max(residPlot);
                    residPlot2 = residPlot;
                    if length(MRS_struct.p.target) == 3 && all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
                        residPlot2(weightRange) = NaN;
                        hold on;
                        plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b' , ...
                            freq(freqbounds), GSHgaussModel(GaussModelParam,freq(freqbounds)), 'r', ...
                            freq(freqbounds), residPlot2, 'k');
                        plot(freq(freqbounds(ChoRange)), residPlot(ChoRange), 'Color', [255 160 64]/255);
                        hold off;
                    else
                        plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b' , ...
                            freq(freqbounds), GSHgaussModel(GaussModelParam,freq(freqbounds)), 'r', ...
                            freq(freqbounds), residPlot2, 'k');
                    end
                    set(gca,'XLim',[1.8 4.2]);
                    
                case 'Lac'
                    plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b', ...
                        freq(freqbounds), FourGaussModel(FourGaussModelParam,freq(freqbounds)), 'r', ...
                        freq(freqbounds), FourGaussModel(MMGaussModelParam,freq(freqbounds)), 'r' , ...
                        freq(freqbounds), resid, 'k');
                    set(gca,'XLim',[0 2.1]);
                    
                case 'Glx'
                    plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b', ...
                        freq(freqbounds), DoubleGaussModel(GaussModelParam,freq(freqbounds)), 'r', ...
                        freq(freqbounds), resid, 'k');
                    set(gca,'XLim',[3.4 4.2]);
                    
                case 'GABAGlx'
                    residPlot = residPlot + metabmin - max(residPlot);
                    residPlot2 = residPlot;
                    residPlot2(weightRange) = NaN;
                    hold on;
                    plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b', ...
                        freq(freqbounds), GABAGlxModel(GaussModelParam,freq(freqbounds)), 'r', ...
                        freq(freqbounds), residPlot2, 'k');
                    % Plot weighted portion of residuals in different color
                    if MRS_struct.p.HERMES && any(strcmp(MRS_struct.p.vendor,{'Philips','Philips_data','Philips_raw'}))
                        plot(freq(freqbounds(ChoRange)), residPlot(ChoRange), 'Color', [255 160 64]/255);
                        plot(freq(freqbounds(GlxDownFieldRange)), residPlot(GlxDownFieldRange), 'Color', [255 160 64]/255);
                    else
                        plot(freq(freqbounds(ChoRange)), residPlot(ChoRange), 'Color', [255 160 64]/255);
                    end
                    hold off;
                    set(gca,'XLim',[2.7 4.2]);
                    
                case 'EtOH'
                    residPlot = residPlot + metabmin - max(residPlot);
                    residPlot2 = residPlot;
                    residPlot2(weightRange) = NaN;
                    hold on;
                    plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b', ...
                        freq(freqbounds), EtOHModel(LorentzModelParam,freq(freqbounds)), 'r', ...
                        freq(freqbounds), residPlot2, 'k');
                    plot(freq(freqbounds(LacRange)), residPlot(LacRange), 'Color', [255 160 64]/255);
                    hold off;
                    set(gca,'XLim',[0.4 1.9]);
            end
            
            % From here on is cosmetic - adding labels etc.
            switch target{jj}
                case 'GABA'
                    text(3, metabmax/4, target{jj}, 'HorizontalAlignment', 'center');
                    labelbounds = freq <= 2.4 & freq >= 2;
                    tailtop = max(real(DIFF(ii,labelbounds)));
                    tailbottom = min(real(DIFF(ii,labelbounds)));
                    text(2.8, min(resid), 'residual', 'HorizontalAlignment', 'left');
                    text(2.8, tailtop+metabmax/20, 'data', 'Color', [0 0 1]);
                    text(2.8, tailbottom-metabmax/20, 'model', 'Color', [1 0 0]);
                    if length(MRS_struct.p.target) ~= 3
                        text(3.2, min(residPlot)-0.5*abs(max(residPlot)), 'weighted', 'Color', [255 160 64]/255, 'HorizontalAlignment', 'center');
                    end
                    
                case 'GSH'
                    text(2.95, maxinGSH+maxinGSH/4, target{jj}, 'HorizontalAlignment', 'center');
                    labelbounds = freq <= 2.4 & freq >= 1.75;
                    tailtop = max(real(DIFF(ii,labelbounds)));
                    tailbottom = min(real(DIFF(ii,labelbounds)));
                    text(2.25, min(resid),'residual', 'HorizontalAlignment', 'left');
                    text(2.25, tailtop+metabmax/20, 'data', 'Color', [0 0 1]);
                    text(2.45, tailbottom-20*metabmax/20, 'model', 'Color', [1 0 0]);
                    if length(MRS_struct.p.target) == 3 && all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
                        text(3.3, max(residPlot)+0.5*abs(max(residPlot)-min(residPlot)), 'weighted', 'Color', [255 160 64]/255, 'HorizontalAlignment', 'right');
                    end
                    
                case 'Glx'
                    text(3.8, metabmax/4, target{jj}, 'HorizontalAlignment', 'center');
                    labelbounds = freq <= 3.6 & freq >= 3.4;
                    tailtop = max(real(DIFF(ii,labelbounds)));
                    tailbottom = min(real(DIFF(ii,labelbounds)));
                    text(3.5, min(resid),'residual', 'HorizontalAlignment', 'left');
                    text(3.5, tailtop+metabmax/20, 'data', 'Color', [0 0 1]);
                    text(3.5, tailbottom-metabmax/20, 'model', 'Color', [1 0 0]);
                    
                case 'GABAGlx'
                    text(3, metabmax/4, 'GABA', 'HorizontalAlignment', 'center');
                    text(3.755, metabmax/4, 'Glx', 'HorizontalAlignment', 'center');
                    text(2.8, min(resid), 'residual', 'HorizontalAlignment', 'left');
                    labelbounds = freq <= 2.8 & freq >= 2.7;
                    tailtop = max(real(DIFF(ii,labelbounds)));
                    tailbottom = min(real(DIFF(ii,labelbounds)));
                    text(2.8, tailtop+metabmax/20, 'data', 'Color', [0 0 1]);
                    text(2.8, tailbottom-metabmax/20, 'model', 'Color', [1 0 0]);
                    text(3.2, min(residPlot)-0.5*abs(max(residPlot)), 'weighted', 'Color', [255 160 64]/255, 'HorizontalAlignment', 'center');
                    
                case 'EtOH'
                    text(1.45,metabmax/4,target{jj}, 'HorizontalAlignment', 'center');
                    labelbounds = freq <= 0.9 & freq >= 0.5;
                    tailtop = max(real(DIFF(ii,labelbounds)));
                    tailbottom = min(real(DIFF(ii,labelbounds)));
                    text(0.6, min(resid),'residual', 'HorizontalAlignment', 'left');
                    text(0.8, tailtop+metabmax/20, 'data', 'Color', [0 0 1]);
                    text(0.8, tailbottom-metabmax/20, 'model', 'Color', [1 0 0]);
                    text(1.4, max(residPlot2)+0.5*abs(max(residPlot2)), 'weighted', 'Color', [255 160 64]/255, 'HorizontalAlignment', 'center');
            end
            
            title('Edited Spectrum and Model Fit');
            xlabel('ppm');
            set(gca,'XDir','reverse','TickDir','out','Box','off');
            ax = get(gca,'YAxis');
            set(ax,'Visible','off');
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   2.  Water Fit
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if strcmp(MRS_struct.p.Reference_compound,'H2O')
                
                % Estimate height and baseline from data
                [maxinWater, watermaxindex] = max(real(WaterData(ii,:)),[],2);
                waterbase = mean(real(WaterData(ii,freq <= 4 & freq >= 3.8)));
                
                % Philips data do not phase well based on first point, so do a preliminary
                % fit, then adjust phase of WaterData accordingly
                % Do this only if the eddy current correction is not performed on water -- MGSaleh 2016
                if strcmpi(MRS_struct.p.vendor,'Philips') && ~MRS_struct.p.water_phase_correction
                    %Run preliminary fit of data
                    LGModelInit = [maxinWater 20 freq(watermaxindex) 0 waterbase -50];
                    %Fit from 5.6 ppm to 3.8 ppm RE 110826
                    freqbounds = freq <= 5.6 & freq >= 3.8;
                    
                    % Do the water fit (Lorentz-Gauss)
                    LGModelParam = nlinfit(freq(freqbounds), real(WaterData(ii,freqbounds)), @LorentzGaussModel, LGModelInit, nlinopts);
                    
                    Eerror = zeros([120 1]);
                    for ll = 1:120
                        Data = WaterData(ii,freqbounds)*exp(1i*pi/180*ll*3);
                        Model = LorentzGaussModel(LGModelParam,freq(freqbounds));
                        Eerror(ll) = sum((real(Data)-Model).^2);
                    end
                    [~,index] = min(Eerror);
                    WaterData(ii,:) = WaterData(ii,:) * exp(1i*pi/180*index*3);
                end
                
                LGPModelInit = [maxinWater 20 freq(watermaxindex) 0 waterbase -50 0];
                lb = [0.01*maxinWater 1 4.6 0 0 -50 -pi];
                ub = [40*maxinWater 100 4.8 0.000001 1 0 pi];
                                
                freqbounds = freq <= 5.6 & freq >= 3.8;
                
                % Least-squares model fitting
                LGPModelInit = lsqcurvefit(@LorentzGaussModelP, LGPModelInit, freq(freqbounds), real(WaterData(ii,freqbounds)), lb, ub, lsqopts);
                [LGPModelParam, residw] = nlinfit(freq(freqbounds), real(WaterData(ii,freqbounds)), @LorentzGaussModelP, LGPModelInit, nlinopts);
                
                WaterArea = sum(real(LorentzGaussModel(LGPModelParam(1:end-1),freq(freqbounds))) - BaselineModel(LGPModelParam(3:5),freq(freqbounds)),2);
                MRS_struct.out.(vox{kk}).water.Area(ii) = WaterArea * abs(freq(1)-freq(2));
                waterheight = LGPModelParam(1);
                MRS_struct.out.(vox{kk}).water.FitError(ii) = 100*std(residw)/waterheight;
                
                LG = real(LorentzGaussModel(LGPModelParam(1:end-1),freq(freqbounds))) - BaselineModel(LGPModelParam(3:5),freq(freqbounds));
                LG = LG./max(LG);
                ind = find(LG >= 0.5);
                f = freq(freqbounds);
                w = abs(f(ind(1)) - f(ind(end)));
                MRS_struct.out.(vox{kk}).water.FWHM(ii) = w * MRS_struct.p.LarmorFreq(ii);
                MRS_struct.out.(vox{kk}).water.ModelParam(ii,:) = LGPModelParam;
                MRS_struct.out.(vox{kk}).water.Resid(ii,:) = residw;
                
                % Calculate SNR of water signal
                noiseSigma_Water = CalcNoise(freq, WaterData(ii,:));
                MRS_struct.out.(vox{kk}).water.SNR(ii) = abs(waterheight)/noiseSigma_Water;
                
                % Water spectrum plot
                hb = subplot(2,2,3);
                watmin = min(real(WaterData(ii,:)));
                watmax = max(real(WaterData(ii,:)));
                residw = residw + watmin - max(residw);
                plot(freq(freqbounds), real(WaterData(ii,freqbounds)), 'b', ...
                    freq(freqbounds), real(LorentzGaussModelP(LGPModelParam,freq(freqbounds))), 'r', ...
                    freq(freqbounds), residw, 'k');
                set(gca,'XDir','reverse','TickDir','out','Box','off','XTick',4.2:0.2:5.2);
                xlim([4.2 5.2]);
                ax = get(gca,'YAxis');
                set(ax,'Visible','off');
                % Add on some labels
                text(4.8, watmax/2, 'Water', 'HorizontalAlignment', 'right');
                labelfreq = freq(freqbounds);
                rlabelbounds = labelfreq <= 4.4 & labelfreq >= 4.25;
                axis_bottom = axis;
                text(4.4, max(min(residw(rlabelbounds))-0.05*watmax, axis_bottom(3)), 'residual', 'HorizontalAlignment', 'left');
                
                % Root sum square fit error and concentration in institutional units
                switch target{jj}
                    case 'GABA'
                        MRS_struct.out.(vox{kk}).GABA.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).GABA.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                        MRS_struct = CalcIU(MRS_struct, vox{kk}, 'GABA', ii);
                        
                    case 'Glx'
                        MRS_struct.out.(vox{kk}).Glx.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).Glx.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                        MRS_struct = CalcIU(MRS_struct, vox{kk}, 'Glx', ii);
                        
                    case 'GABAGlx'
                        MRS_struct.out.(vox{kk}).GABA.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).GABA.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                        MRS_struct.out.(vox{kk}).Glx.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).Glx.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                        MRS_struct = CalcIU(MRS_struct, vox{kk}, 'GABA', ii);
                        MRS_struct = CalcIU(MRS_struct, vox{kk}, 'Glx', ii);
                        
                    case 'GSH'
                        MRS_struct.out.(vox{kk}).GSH.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).GSH.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                        MRS_struct = CalcIU(MRS_struct, vox{kk}, (target{jj}), ii);
                        
                    case 'Lac'
                        MRS_struct.out.(vox{kk}).Lac.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).Lac.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                        MRS_struct = CalcIU(MRS_struct, vox{kk}, (target{jj}), ii);
                        
                    case 'EtOH'
                        MRS_struct.out.(vox{kk}).EtOH.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).EtOH.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                        MRS_struct = CalcIU(MRS_struct, vox{kk}, (target{jj}), ii);
                end
                
                % Generate scaled spectra (for plotting)
                MRS_struct.spec.(vox{kk}).(target{jj}).off_scaled(ii,:) = ...
                    MRS_struct.spec.(vox{kk}).(target{jj}).off(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
                if MRS_struct.p.HERMES
                    MRS_struct.spec.(vox{kk}).(target{jj}).off_off_scaled(ii,:) = ...
                        MRS_struct.spec.(vox{kk}).(target{jj}).off_off(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
                end
                MRS_struct.spec.(vox{kk}).(target{jj}).on_scaled(ii,:) = ...
                    MRS_struct.spec.(vox{kk}).(target{jj}).on(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
                MRS_struct.spec.(vox{kk}).(target{jj}).diff_scaled(ii,:) = ...
                    MRS_struct.spec.(vox{kk}).(target{jj}).diff(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
                
                % Reorder structure fields
                MRS_struct.out.(vox{kk}).water = orderfields(MRS_struct.out.(vox{kk}).water, {'Area', 'FWHM', 'SNR', 'ModelParam', 'Resid', 'FitError'});
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   2a.  Residual Water Fit
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if MRS_struct.p.FitResidWater && ~MRS_struct.p.HERMES
                    
                    water_OFF = OFF(ii,:);
                    freqWaterOFF = freq <= 4.68+0.4 & freq >= 4.68-0.4;
                    water_OFF = water_OFF(freqWaterOFF);
                    freqWaterOFF = freq(freqWaterOFF);
                    
                    [maxResidWater, maxInd] = max(abs(real(water_OFF)));
                    s = sign(real(water_OFF(maxInd)));
                    maxResidWater = s * maxResidWater;
                    offset = real(water_OFF(1));
                    
                    LGPModelInit = [maxResidWater 25 freqWaterOFF(maxInd) 0 offset 0.001 0];
                    LGPModelInit([1 5]) = LGPModelInit([1 5]) / maxResidWater; % Scale initial conditions to avoid warnings about numerical underflow
                    
                    %lb = [maxResidWater-abs(2*maxResidWater) 1 freqWaterOFF(maxInd)-0.2 0 0 -200 -pi];
                    %ub = [maxResidWater+abs(2*maxResidWater) 100 freqWaterOFF(maxInd)+0.2 0.000001 1 0 pi];
                    %lb([1 4 5]) = lb([1 4 5]) / maxResidWater;
                    %ub([1 4 5]) = ub([1 4 5]) / maxResidWater;
                    
                    % Least-squares model fitting
                    %LGPModelInit = lsqcurvefit(@LorentzGaussModelP, LGPModelInit, freqWaterOFF, real(water_OFF) / maxResidWater, lb, ub, lsqopts);
                    [MRS_struct.out.(vox{kk}).ResidWater.ModelParam(ii,:), residRW] = nlinfit(freqWaterOFF, real(water_OFF) / maxResidWater, @LorentzGaussModelP, LGPModelInit, nlinopts);
                    
                    % Rescale fit parameters and residuals
                    MRS_struct.out.(vox{kk}).ResidWater.ModelParam(ii,[1 4 5]) = MRS_struct.out.(vox{kk}).ResidWater.ModelParam(ii,[1 4 5]) * maxResidWater;
                    residRW = residRW * maxResidWater;
                    
                    MRS_struct.out.(vox{kk}).ResidWater.FitError(ii) = 100*std(residRW)/MRS_struct.out.(vox{kk}).ResidWater.ModelParam(ii,1);
                    
                    MRS_struct.out.(vox{kk}).ResidWater.SuppressionFactor(ii) = ...
                        (MRS_struct.out.(vox{kk}).water.ModelParam(ii,1) - abs(MRS_struct.out.(vox{kk}).ResidWater.ModelParam(ii,1))) ...
                        / MRS_struct.out.(vox{kk}).water.ModelParam(ii,1);
                    
                end
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   3.  Cr Fit
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Cr_OFF = OFF(ii,:);
            freqboundsChoCr = freq <= 3.6 & freq >= 2.6;
            
            % Do some detective work to figure out the initial parameters
            ChoCrMeanSpec = Cr_OFF(freqboundsChoCr).';
            Baseline_offset = real(ChoCrMeanSpec(1)+ChoCrMeanSpec(end))/2;
            Width_estimate = 0.05;
            Area_estimate = (max(real(ChoCrMeanSpec))-min(real(ChoCrMeanSpec)))*Width_estimate*4;
            ChoCr_initx = [Area_estimate Width_estimate 3.02 0 Baseline_offset 0 1] ...
                .* [1 2*MRS_struct.p.LarmorFreq(ii) MRS_struct.p.LarmorFreq(ii) 180/pi 1 1 1];
            ChoCrModelParam = FitChoCr(freq(freqboundsChoCr), ChoCrMeanSpec, ChoCr_initx, MRS_struct.p.LarmorFreq(ii));
            MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,:) = ChoCrModelParam ./ [1 2*MRS_struct.p.LarmorFreq(ii) MRS_struct.p.LarmorFreq(ii) 180/pi 1 1 1];
            
            % Initialise fitting pars
            freqboundsCr = freq <= 3.12 & freq >= 2.72;
            LorentzModelInit = [max(real(Cr_OFF(freqboundsCr))) 0.05 3.0 0 0 0];
            
            % Least-squares model fitting
            LorentzModelInit = lsqcurvefit(@LorentzModel, LorentzModelInit, freq(freqboundsCr), real(Cr_OFF(freqboundsCr)), [], [], lsqopts);
            [LorentzModelParam, residCr] = nlinfit(freq(freqboundsCr), real(Cr_OFF(freqboundsCr)), @LorentzModel, LorentzModelInit, nlinopts);
            
            MRS_struct.out.(vox{kk}).Cr.ModelParam(ii,:) = LorentzModelParam;
            Crheight = LorentzModelParam(1)/(2*pi*LorentzModelParam(2));
            MRS_struct.out.(vox{kk}).Cr.FitError(ii) = 100*std(residCr)/Crheight;
            MRS_struct.out.(vox{kk}).Cr.Resid(ii,:) = residCr;
            
            MRS_struct.out.(vox{kk}).Cr.Area(ii) = sum(real(TwoLorentzModel([MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,1:end-1) 0],freq(freqboundsChoCr)) - ...
                TwoLorentzModel([0 MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,2:end-1) 0],freq(freqboundsChoCr)))) * abs(freq(1)-freq(2));
            MRS_struct.out.(vox{kk}).Cho.Area(ii) = sum(real(TwoLorentzModel(MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,:),freq(freqboundsChoCr)) - ...
                TwoLorentzModel([MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,1:end-1) 0],freq(freqboundsChoCr)))) * abs(freq(1)-freq(2));
            MRS_struct.out.(vox{kk}).Cr.FWHM(ii) = ChoCrModelParam(2);
            
            % Calculate SNR of Cr signal
            noiseSigma_OFF = CalcNoise(freq, OFF(ii,:));
            MRS_struct.out.(vox{kk}).Cr.SNR(ii) = abs(Crheight)/noiseSigma_OFF;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   4.  NAA Fit
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            NAA_OFF = OFF(ii,:);
            freqbounds = find(freq <= 2.25 & freq >= 1.75);
            
            maxinNAA = max(real(NAA_OFF(freqbounds)));
            grad_points = (real(NAA_OFF(freqbounds(end))) - real(NAA_OFF(freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1)); %in points
            LinearInit = grad_points ./ abs(freq(1) - freq(2));
            constInit = (real(NAA_OFF(freqbounds(end))) + real(NAA_OFF(freqbounds(1)))) ./ 2;
            
            LorentzModelInit = [maxinNAA 0.05 2.01 0 -LinearInit constInit];
            lb = [0 0.01 1.97 0 -40*maxinNAA -2000*maxinNAA];
            ub = [4000*maxinNAA 0.1 2.05 0.5 40*maxinNAA 1000*maxinNAA];
            
            % Least-squares model fitting
            LorentzModelInit = lsqcurvefit(@LorentzModel, LorentzModelInit, freq(freqbounds), real(NAA_OFF(freqbounds)), lb, ub, lsqopts);
            [LorentzModelParam, resid] = nlinfit(freq(freqbounds), real(NAA_OFF(freqbounds)), @LorentzModel, LorentzModelInit, nlinopts);
            
            NAAheight = LorentzModelParam(1)/(2*pi*LorentzModelParam(2));
            MRS_struct.out.(vox{kk}).NAA.FitError(ii) = 100*std(resid)/NAAheight;
            NAAModelParam = LorentzModelParam;
            NAAModelParam(4) = 0;
            MRS_struct.out.(vox{kk}).NAA.Area(ii) = sum(LorentzModel(NAAModelParam,freq(freqbounds)) - BaselineModel(NAAModelParam([3 6 5]),freq(freqbounds)), 2) * abs(freq(1)-freq(2));
            MRS_struct.out.(vox{kk}).NAA.FWHM(ii) = abs((2*MRS_struct.p.LarmorFreq(ii)) * NAAModelParam(2));
            MRS_struct.out.(vox{kk}).NAA.ModelParam(ii,:) = LorentzModelParam;
            MRS_struct.out.(vox{kk}).NAA.Resid(ii,:) = resid;
            
            % Calculate SNR of NAA signal
            MRS_struct.out.(vox{kk}).NAA.SNR(ii) = abs(NAAheight) / noiseSigma_OFF;
            
            % Root sum square fit errors and concentrations as metabolite ratios
            if strcmpi(target{jj},'GABAGlx')
                MRS_struct.out.(vox{kk}).GABA.FitError_Cr(ii) = sqrt(MRS_struct.out.(vox{kk}).GABA.FitError(ii).^2 + MRS_struct.out.(vox{kk}).Cr.FitError(ii).^2);
                MRS_struct.out.(vox{kk}).Glx.FitError_Cr(ii) = sqrt(MRS_struct.out.(vox{kk}).Glx.FitError(ii).^2 + MRS_struct.out.(vox{kk}).Cr.FitError(ii).^2);
                MRS_struct.out.(vox{kk}).GABA.FitError_NAA(ii) = sqrt(MRS_struct.out.(vox{kk}).GABA.FitError(ii).^2 + MRS_struct.out.(vox{kk}).NAA.FitError(ii).^2);
                MRS_struct.out.(vox{kk}).Glx.FitError_NAA(ii) = sqrt(MRS_struct.out.(vox{kk}).Glx.FitError(ii).^2 + MRS_struct.out.(vox{kk}).NAA.FitError(ii).^2);
                MRS_struct.out.(vox{kk}).GABA.ConcCr(ii) = MRS_struct.out.(vox{kk}).GABA.Area(ii) / MRS_struct.out.(vox{kk}).Cr.Area(ii);
                MRS_struct.out.(vox{kk}).GABA.ConcCho(ii) = MRS_struct.out.(vox{kk}).GABA.Area(ii) / MRS_struct.out.(vox{kk}).Cho.Area(ii);
                MRS_struct.out.(vox{kk}).GABA.ConcNAA(ii) = MRS_struct.out.(vox{kk}).GABA.Area(ii) / MRS_struct.out.(vox{kk}).NAA.Area(ii);
                MRS_struct.out.(vox{kk}).Glx.ConcCr(ii) = MRS_struct.out.(vox{kk}).Glx.Area(ii) / MRS_struct.out.(vox{kk}).Cr.Area(ii);
                MRS_struct.out.(vox{kk}).Glx.ConcCho(ii) = MRS_struct.out.(vox{kk}).Glx.Area(ii) / MRS_struct.out.(vox{kk}).Cho.Area(ii);
                MRS_struct.out.(vox{kk}).Glx.ConcNAA(ii) = MRS_struct.out.(vox{kk}).Glx.Area(ii) / MRS_struct.out.(vox{kk}).NAA.Area(ii);
            else
                MRS_struct.out.(vox{kk}).(target{jj}).FitError_Cr(ii) = sqrt(MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii).^2 + MRS_struct.out.(vox{kk}).Cr.FitError(ii).^2);
                MRS_struct.out.(vox{kk}).(target{jj}).FitError_NAA(ii) = sqrt(MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii).^2 + MRS_struct.out.(vox{kk}).NAA.FitError(ii).^2);
                MRS_struct.out.(vox{kk}).(target{jj}).ConcCr(ii) = MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) / MRS_struct.out.(vox{kk}).Cr.Area(ii);
                MRS_struct.out.(vox{kk}).(target{jj}).ConcCho(ii) = MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) / MRS_struct.out.(vox{kk}).Cho.Area(ii);
                MRS_struct.out.(vox{kk}).(target{jj}).ConcNAA(ii) = MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) / MRS_struct.out.(vox{kk}).NAA.Area(ii);
            end
            
            % Reorder structure fields
            if ~MRS_struct.p.HERMES
                if strcmp(MRS_struct.p.Reference_compound,'H2O')
                    if strcmpi(target{jj},'GABAGlx')
                        MRS_struct.out.(vox{kk}).GABA = orderfields(MRS_struct.out.(vox{kk}).GABA, ...
                            {'Area', 'FWHM', 'SNR', 'ModelParam', 'Resid', 'FitError', 'FitError_W', 'FitError_Cr', 'FitError_NAA', 'ConcIU', 'ConcCr', 'ConcCho', 'ConcNAA'});
                        MRS_struct.out.(vox{kk}).Glx = orderfields(MRS_struct.out.(vox{kk}).Glx, ...
                            {'Area', 'FWHM', 'SNR', 'ModelParam', 'Resid', 'FitError', 'FitError_W', 'FitError_Cr', 'FitError_NAA', 'ConcIU', 'ConcCr', 'ConcCho', 'ConcNAA'});
                    else
                        MRS_struct.out.(vox{kk}).(target{jj}) = orderfields(MRS_struct.out.(vox{kk}).(target{jj}), ...
                            {'Area', 'FWHM', 'SNR', 'ModelParam', 'Resid', 'FitError', 'FitError_W', 'FitError_Cr', 'FitError_NAA', 'ConcIU', 'ConcCr', 'ConcCho', 'ConcNAA'});
                    end
                else
                    if strcmpi(target{jj},'GABAGlx')
                        MRS_struct.out.(vox{kk}).GABA = orderfields(MRS_struct.out.(vox{kk}).GABA, ...
                            {'Area', 'FWHM', 'SNR', 'ModelParam', 'Resid', 'FitError', 'FitError_Cr', 'FitError_NAA', 'ConcCr', 'ConcCho', 'ConcNAA'});
                        MRS_struct.out.(vox{kk}).Glx = orderfields(MRS_struct.out.(vox{kk}).Glx, ...
                            {'Area', 'FWHM', 'SNR', 'ModelParam', 'Resid', 'FitError', 'FitError_Cr', 'FitError_NAA', 'ConcCr', 'ConcCho', 'ConcNAA'});
                    else
                        MRS_struct.out.(vox{kk}).(target{jj}) = orderfields(MRS_struct.out.(vox{kk}).(target{jj}), ...
                            {'Area', 'FWHM', 'SNR', 'ModelParam', 'Resid', 'FitError', 'FitError_Cr', 'FitError_NAA', 'ConcCr', 'ConcCho', 'ConcNAA'});
                    end
                end
            end
            MRS_struct.out.(vox{kk}).Cr = orderfields(MRS_struct.out.(vox{kk}).Cr, {'Area', 'FWHM', 'SNR', 'ModelParam', 'Resid', 'FitError'});
            MRS_struct.out.(vox{kk}).NAA = orderfields(MRS_struct.out.(vox{kk}).NAA, {'Area', 'FWHM', 'SNR', 'ModelParam', 'Resid', 'FitError'});
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   5. Build GannetFit Output
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Crmin = min(real(Cr_OFF(freqboundsCr)));
            Crmax = max(real(Cr_OFF(freqboundsCr)));
            resmaxCr = max(residCr);
            residCr = residCr + Crmin - resmaxCr;
            if strcmp(MRS_struct.p.Reference_compound,'H2O')
                hd = subplot(2,2,4);
                plot(freq, real(Cr_OFF), 'b', ...
                    freq(freqboundsChoCr), real(TwoLorentzModel(MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,:),freq(freqboundsChoCr))), 'r', ...
                    freq(freqboundsChoCr), real(TwoLorentzModel([MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,1:(end-1)) 0],freq(freqboundsChoCr))), 'r', ...
                    freq(freqboundsCr), residCr, 'k');
                set(gca,'XDir','reverse','TickDir','out','Box','off','XTick',2.6:0.2:3.6);
                xlim([2.6 3.6]);
                ax = get(gca,'YAxis');
                set(ax,'Visible','off');
                xlabel('ppm');
                text(2.94, Crmax*0.75, 'Creatine', 'HorizontalAlignment', 'left');
                subplot(2,2,3);
                [~,hi] = inset(hb,hd);
                set(hi,'fontsize',6);
                insert = get(hi,'pos');
                axi = get(hb,'pos');
                set(hi,'pos',[axi(1)+axi(3)-insert(3) insert(2:4)]);
                text(4.8, watmax/2, 'Water', 'HorizontalAlignment', 'right');
                set(gca,'XDir','reverse','TickDir','out','Box','off');
                ax = get(gca,'YAxis');
                set(ax,'Visible','off');
                xlabel('ppm');
                title('Reference Signals');
            else
                hb = subplot(2,2,3);
                plot(freq, real(Cr_OFF), 'b', ...
                    freq(freqboundsChoCr), real(TwoLorentzModel(MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,:),freq(freqboundsChoCr))), 'r', ...
                    freq(freqboundsChoCr), real(TwoLorentzModel([MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,1:(end-1)) 0],freq(freqboundsChoCr))), 'r', ...
                    freq(freqboundsCr), residCr, 'k');
                set(gca,'XDir','reverse','TickDir','out','Box','off','XTick',2.6:0.2:3.6);
                ax = get(gca,'YAxis');
                set(ax,'Visible','off');
                xlim([2.6 3.6]);
                crlabelbounds = freq(freqboundsCr) <= 3.12 & freq(freqboundsCr) >= 2.72;
                hcres = text(3, max(residCr(crlabelbounds))+0.05*Crmax, 'residual');
                set(hcres,'HorizontalAlignment', 'center');
                text(2.7,0.1*Crmax,'data','Color',[0 0 1]);
                text(2.7,0.01*Crmax,'model','Color',[1 0 0]);
                text(2.94,Crmax*0.75,'Creatine');
                xlabel('ppm');
                title('Reference Signal');
            end
            
            if any(strcmp('mask',fieldnames(MRS_struct)))
                hc = subplot(2,2,2);
                set(hc,'pos',[0.52 0.52 0.42 0.42]) % move the axes slightly
                size_max = size(MRS_struct.mask.img{ii},1);
                imagesc(MRS_struct.mask.img{ii}(:,size_max+(1:size_max)));
                colormap('gray');
                caxis([0 1])
                axis equal;
                axis tight;
                axis off;
                subplot(2,2,4,'replace');
            else
                subplot(3,2,[2 4]);
                axis off;
            end
            
            % Cleaner text alignment; move GABA/Glx to separate lines
            text_pos = 0.95; % A variable to determine y-position of text on printout on figure
            shift = 0.06;
            
            % 1. Filename
            if strcmp(MRS_struct.p.vendor,'Siemens_rda')
                [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii*2-1});
            else
                [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii});
            end
            fname = [tmp tmp2];
            if length(fname) > 30
                fname = [fname(1:30) '...'];
            end
            
            text(0.4, text_pos, 'Filename: ', 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
            text(0.425, text_pos, fname, 'FontName', 'Arial', 'FontSize', 10, 'Interpreter', 'none');
            
            % 2a. Area
            text(0.4, text_pos-shift, 'Area  ', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
            
            switch target{jj}
                case 'GABA'
                    tmp2 = 'GABA+: ';
                    tmp3 = sprintf('%.3g', MRS_struct.out.(vox{kk}).GABA.Area(ii));
                case 'Glx'
                    tmp2 = 'Glx: ';
                    tmp3 = sprintf('%.3g', MRS_struct.out.(vox{kk}).Glx.Area(ii));
                case 'GABAGlx'
                    tmp2 = 'GABA+: ';
                    tmp3 = sprintf('%.3g', MRS_struct.out.(vox{kk}).GABA.Area(ii));
                    tmp4 = 'Glx: ';
                    tmp5 = sprintf('%.3g', MRS_struct.out.(vox{kk}).Glx.Area(ii));
                case 'GSH'
                    tmp2 = 'GSH: ';
                    tmp3 = sprintf('%.3g', MRS_struct.out.(vox{kk}).(target{jj}).Area(ii));
                case 'Lac'
                    tmp2 = 'Lac: ';
                    tmp3 = sprintf('%.3g', MRS_struct.out.(vox{kk}).(target{jj}).Area(ii));
                case 'EtOH'
                    tmp2 = 'EtOH: ';
                    tmp3 = sprintf('%.3g', MRS_struct.out.(vox{kk}).(target{jj}).Area(ii));
            end
            if strcmp(target{jj},'GABAGlx')
                text(0.4, text_pos-2*shift, tmp2, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                text(0.425, text_pos-2*shift, tmp3, 'FontName', 'Arial', 'FontSize', 10);
                text(0.4, text_pos-3*shift, tmp4, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                text(0.425, text_pos-3*shift, tmp5, 'FontName', 'Arial', 'FontSize', 10);
                n = 0;
            else
                text(0.4, text_pos-2*shift, tmp2, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                text(0.425, text_pos-2*shift, tmp3, 'FontName', 'Arial', 'FontSize', 10);
                n = shift;
            end
            
            if strcmp(MRS_struct.p.Reference_compound,'H2O')
                
                % 2b. Area (Water / Cr)
                tmp1 = 'Water: ';
                tmp2 = sprintf('%.3g', MRS_struct.out.(vox{kk}).water.Area(ii));
                tmp3 = 'Cr: ';
                tmp4 = sprintf('%.3g', MRS_struct.out.(vox{kk}).Cr.Area(ii));
                
                text(0.4, text_pos-4*shift+n, tmp1, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                text(0.425, text_pos-4*shift+n, tmp2, 'FontName', 'Arial', 'FontSize', 10);
                text(0.4, text_pos-5*shift+n, tmp3, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                text(0.425, text_pos-5*shift+n, tmp4, 'FontName', 'Arial', 'FontSize', 10);
                
                % 3. FWHM
                tmp1 = 'FWHM  ';
                text(0.4, text_pos-6*shift+n, tmp1, 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                
                tmp2 = sprintf('%.2f Hz', MRS_struct.out.(vox{kk}).water.FWHM(ii));
                tmp3 = sprintf('%.2f Hz', MRS_struct.out.(vox{kk}).Cr.FWHM(ii));
                text(0.4, text_pos-7*shift+n, 'Water: ', 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                text(0.425, text_pos-7*shift+n, tmp2, 'FontName', 'Arial', 'FontSize', 10);
                text(0.4, text_pos-8*shift+n, 'Cr: ', 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                text(0.425, text_pos-8*shift+n, tmp3, 'FontName', 'Arial', 'FontSize', 10);
                
                % 4a. Fit Error
                tmp1 = 'Fit Error  ';
                text(0.4, text_pos-9*shift+n, tmp1, 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                
                switch target{jj}
                    
                    case {'GABA','Glx','GSH','Lac','EtOH'}
                        % 4b. Fit Error
                        if strcmpi(target{jj},'GABA')
                            tmp2 = 'GABA+,Water: ';
                            tmp3 = 'GABA+,Cr: ';
                        else
                            tmp2 = sprintf('%s,Water: ', target{jj});
                            tmp3 = sprintf('%s,Cr: ', target{jj});
                        end
                        tmp4 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).(target{jj}).FitError_W(ii));
                        tmp5 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).(target{jj}).FitError_Cr(ii));
                        
                        text(0.4, text_pos-10*shift+n, tmp2, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-10*shift+n, tmp4, 'FontName', 'Arial', 'FontSize', 10);
                        text(0.4, text_pos-11*shift+n, tmp3, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-11*shift+n, tmp5, 'FontName', 'Arial', 'FontSize', 10);
                        
                        % 5. Quantification
                        tmp1 = 'Quantification';
                        if strcmpi(target{jj},'GABA')
                            tmp2 = 'GABA+/Water: ';
                            tmp4 = 'GABA+/Cr: ';
                        else
                            tmp2 = sprintf('%s/Water: ', target{jj});
                            tmp4 = sprintf('%s/Cr: ', target{jj});
                        end
                        tmp3 = sprintf('%.2f i.u.', MRS_struct.out.(vox{kk}).(target{jj}).ConcIU(ii));
                        tmp5 = sprintf('%.2f', MRS_struct.out.(vox{kk}).(target{jj}).ConcCr(ii));
                        
                        text(0.4, text_pos-12*shift+n, tmp1, 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                        text(0.4, text_pos-13*shift+n, tmp2, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-13*shift+n, tmp3, 'FontName', 'Arial', 'FontSize', 10);
                        text(0.4, text_pos-14*shift+n, tmp4, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-14*shift+n, tmp5, 'FontName', 'Arial', 'FontSize', 10);
                        n = 5*shift;
                        
                    case 'GABAGlx'
                        % 4b. Fit Error
                        tmp2 = 'GABA+,Water: ';
                        tmp3 = 'GABA+,Cr: ';
                        tmp4 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).GABA.FitError_W(ii));
                        tmp5 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).GABA.FitError_Cr(ii));
                        tmp6 = 'Glx,Water: ';
                        tmp7 = 'Glx,Cr: ';
                        tmp8 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).Glx.FitError_W(ii));
                        tmp9 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).Glx.FitError_Cr(ii));
                        
                        text(0.4, text_pos-10*shift, tmp2, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-10*shift, tmp4, 'FontName', 'Arial', 'FontSize', 10);
                        text(0.4, text_pos-11*shift, tmp3, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-11*shift, tmp5, 'FontName', 'Arial', 'FontSize', 10);
                        text(0.4, text_pos-12*shift, tmp6, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-12*shift, tmp8, 'FontName', 'Arial', 'FontSize', 10);
                        text(0.4, text_pos-13*shift, tmp7, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-13*shift, tmp9, 'FontName', 'Arial', 'FontSize', 10);
                        
                        % 5. Quantification
                        tmp1 = 'Quantification  ';
                        tmp2 = 'GABA+/Water: ';
                        tmp3 = sprintf('%.2f i.u.', MRS_struct.out.(vox{kk}).GABA.ConcIU(ii));
                        tmp4 = 'GABA+/Cr: ';
                        tmp5 = sprintf('%.2f', MRS_struct.out.(vox{kk}).GABA.ConcCr(ii));
                        tmp6 = 'Glx/Water: ';
                        tmp7 = sprintf('%.2f i.u.', MRS_struct.out.(vox{kk}).Glx.ConcIU(ii));
                        tmp8 = 'Glx/Cr: ';
                        tmp9 = sprintf('%.2f', MRS_struct.out.(vox{kk}).Glx.ConcCr(ii));
                        
                        text(0.4, text_pos-14*shift, tmp1, 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                        text(0.4, text_pos-15*shift, tmp2, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-15*shift, tmp3, 'FontName', 'Arial', 'FontSize', 10);
                        text(0.4, text_pos-16*shift, tmp4, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-16*shift, tmp5, 'FontName', 'Arial', 'FontSize', 10);
                        text(0.4, text_pos-17*shift, tmp6, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-17*shift, tmp7, 'FontName', 'Arial', 'FontSize', 10);
                        text(0.4, text_pos-18*shift, tmp8, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-18*shift, tmp9, 'FontName', 'Arial', 'FontSize', 10);
                        n = 0;
                        
                end
                
                % 6. FitVer
                text(0.4, text_pos-19*shift+n, 'FitVer: ', 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                text(0.425, text_pos-19*shift+n, MRS_struct.version.fit, 'FontName', 'Arial', 'FontSize', 10);
                
            else
                
                tmp1 = 'Cr: ';
                tmp2 = sprintf('%.3g', MRS_struct.out.(vox{kk}).Cr.Area(ii));
                
                text(0.4, text_pos-4*shift+n, tmp1, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                text(0.425, text_pos-4*shift+n, tmp2, 'FontName', 'Arial', 'FontSize', 10);
                
                tmp1 = 'FWHM  ';
                text(0.4, text_pos-5*shift+n, tmp1, 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                
                tmp1 = sprintf('%.2f Hz', MRS_struct.out.(vox{kk}).Cr.FWHM(ii));
                text(0.4, text_pos-6*shift+n, 'Cr: ', 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                text(0.425, text_pos-6*shift+n, tmp1, 'FontName', 'Arial', 'FontSize', 10);
                
                % 4a. Fit Error
                tmp1 = 'Fit Error  ';
                text(0.4, text_pos-7*shift+n, tmp1, 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                
                switch target{jj}
                    
                    case {'GABA','Glx','GSH','Lac','EtOH'}
                        % 4b. Fit Error
                        if strcmpi(target{jj},'GABA')
                            tmp2 = 'GABA+,Cr: ';
                        else
                            tmp2 = sprintf('%s,Cr: ', target{jj});
                        end
                        tmp3 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).(target{jj}).FitError_Cr(ii));
                        
                        text(0.4, text_pos-8*shift+n, tmp2, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-8*shift+n, tmp3, 'FontName', 'Arial', 'FontSize', 10);
                        
                        % 5. Quantification
                        tmp1 = 'Quantification';
                        if strcmpi(target{jj},'GABA')
                            tmp2 = 'GABA+/Cr: ';
                        else
                            tmp2 = sprintf('%s/Cr: ', target{jj});
                        end
                        tmp3 = sprintf('%.2f', MRS_struct.out.(vox{kk}).(target{jj}).ConcCr(ii));
                        
                        text(0.4, text_pos-9*shift+n, tmp1, 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                        text(0.4, text_pos-10*shift+n, tmp2, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-10*shift+n, tmp3, 'FontName', 'Arial', 'FontSize', 10);
                        n = 3*shift;
                        
                    case 'GABAGlx'
                        % 4b. Fit Error
                        tmp1 = 'GABA+,Cr: ';
                        tmp2 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).GABA.FitError_Cr(ii));
                        tmp3 = 'Glx,Cr: ';
                        tmp4 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).Glx.FitError_Cr(ii));
                        
                        text(0.4, text_pos-8*shift, tmp1, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-8*shift, tmp2, 'FontName', 'Arial', 'FontSize', 10);
                        text(0.4, text_pos-9*shift, tmp3, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-9*shift, tmp4, 'FontName', 'Arial', 'FontSize', 10);
                        
                        % 5. Quantification
                        tmp1 = 'Quantification  ';
                        tmp2 = 'GABA+/Cr: ';
                        tmp3 = sprintf('%.2f', MRS_struct.out.(vox{kk}).GABA.ConcCr(ii));
                        tmp4 = 'Glx/Cr: ';
                        tmp5 = sprintf('%.2f', MRS_struct.out.(vox{kk}).Glx.ConcCr(ii));
                        
                        text(0.4, text_pos-10*shift, tmp1, 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                        text(0.4, text_pos-11*shift, tmp2, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-11*shift, tmp3, 'FontName', 'Arial', 'FontSize', 10);
                        text(0.4, text_pos-12*shift, tmp4, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                        text(0.425, text_pos-12*shift, tmp5, 'FontName', 'Arial', 'FontSize', 10);
                        n = 0;
                        
                end
                
                % 6. FitVer
                text(0.4, text_pos-13*shift+n, 'FitVer: ', 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                text(0.425, text_pos-13*shift+n, MRS_struct.version.fit, 'FontName', 'Arial', 'FontSize', 10);
                
            end
            
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
            
            if strcmp(MRS_struct.p.vendor,'Siemens_rda')
                [~,metabfile_nopath] = fileparts(MRS_struct.metabfile{ii*2-1});
            else
                [~,metabfile_nopath] = fileparts(MRS_struct.metabfile{ii});
            end
            
            if any(strcmp(listfonts,'Arial'))
                set(findall(h,'-property','FontName'),'FontName','Arial');
            end
            
            % Create output folder
            if ~exist(fullfile(pwd, 'GannetFit_output'),'dir')
                mkdir(fullfile(pwd, 'GannetFit_output'));
            end
            
            % Save PDF output
            set(h,'PaperUnits','inches');
            set(h,'PaperSize',[11 8.5]);
            set(h,'PaperPosition',[0 0 11 8.5]);
            
            if strcmpi(MRS_struct.p.vendor,'Philips_data')
                pdfname = fullfile(pwd, 'GannetFit_output', [fullpath '_' target{jj} '_' vox{kk} '_fit.pdf']);
            else
                pdfname = fullfile(pwd, 'GannetFit_output', [metabfile_nopath '_' target{jj} '_' vox{kk} '_fit.pdf']);
            end
            saveas(h, pdfname);
            
        end
        
    end
    
    % Reorder structure
    if isfield(MRS_struct, 'mask')
        if isfield(MRS_struct, 'waterfile')
            structorder = {'version', 'ii', 'metabfile', 'waterfile', 'p', 'fids', 'spec', 'out', 'mask'};
        else
            structorder = {'version', 'ii', 'metabfile', 'p', 'fids', 'spec', 'out', 'mask'};
        end
    else
        if isfield(MRS_struct, 'waterfile')
            structorder = {'version', 'ii', 'metabfile', 'waterfile', 'p', 'fids', 'spec', 'out'};
        else
            structorder = {'version','ii', 'metabfile', 'p', 'fids', 'spec', 'out'};
        end
    end
    MRS_struct = orderfields(MRS_struct, structorder);
    
    if MRS_struct.p.mat % save MRS_struct as mat file
        mat_name = ['GannetFit_output/MRS_struct_' vox{kk} '.mat'];
        save(mat_name,'MRS_struct');
    end
    
    if MRS_struct.p.csv % export MRS_struct fields into csv file
        ExportToCSV(MRS_struct, target, kk, 'fit');
    end
    
    fprintf('\n\n');
    
end

warning('on');


%%%%%%%%%%%%%%%% GAUSS MODEL %%%%%%%%%%%%%%%%
function F = GaussModel(x,freq)
% Function for Gauss Model

% x(1) = gaussian amplitude
% x(2) = 1/(2*sigma^2)
% x(3) = centre freq of peak
% x(4) = amplitude of linear baseline
% x(5) = constant amplitude offset

F = x(1)*exp(x(2)*(freq-x(3)).*(freq-x(3))) + x(4)*(freq-x(3)) + x(5);


%%%%%%%%%%%%%%%% LORENTZGAUSSMODEL %%%%%%%%%%%%%%%%
function F = LorentzGaussModel(x,freq)
% Function for LorentzGaussModel Model

% CJE 24Nov10 - removed phase term from fit - this is now dealt with
% by the phasing of the water ref scans in MRSLoadPfiles
% Lorentzian Model multiplied by a Gaussian.
% x(1) = Amplitude of (scaled) Lorentzian
% x(2) = 1 / hwhm of Lorentzian (hwhm = half width at half max)
% x(3) = centre freq of Lorentzian
% x(4) = linear baseline slope
% x(5) = constant baseline amplitude
% x(6) =  -1 / 2 * sigma^2  of gaussian

% Lorentzian  = (1/pi) * (hwhm) / (deltaf^2 + hwhm^2)
% Peak height of Lorentzian = 4 / (pi*hwhm)
% F is a normalised Lorentzian - height independent of hwhm
%   = Lorentzian / Peak

F = (x(1)*ones(size(freq))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1)) ... % Lorentzian
    .* (exp(x(6)*(freq-x(3)).*(freq-x(3)))) ... % Gaussian
    + x(4)*(freq-x(3)) ... % linear baseline
    + x(5); % constant baseline


%%%%%%%%%%%%%%%% LORENTZGAUSSMODEL WITH PHASE %%%%%%%%%%%%%%%%
function F = LorentzGaussModelP(x,freq)
% Function for LorentzGaussModel Model with Phase

% Lorentzian Model multiplied by a Gaussian
% x(1) = Amplitude of (scaled) Lorentzian
% x(2) = 1 / hwhm of Lorentzian (hwhm = half width at half max)
% x(3) = centre freq of Lorentzian
% x(4) = linear baseline slope
% x(5) = constant baseline amplitude
% x(6) =  -1 / 2 * sigma^2  of gaussian
% x(7) = phase (in rad)

% Lorentzian  = (1/pi) * (hwhm) / (deltaf^2 + hwhm^2)
% Peak height of Lorentzian = 4 / (pi*hwhm)
% F is a normalised Lorentzian - height independent of hwhm
%   = Lorentzian / Peak

F = ((cos(x(7))*x(1)*ones(size(freq)) + sin(x(7))*x(1)*x(2)*(freq-x(3)))./(x(2)^2*(freq-x(3)).*(freq-x(3))+1)) ... % Lorentzian
    .* (exp(x(6)*(freq-x(3)).*(freq-x(3)))) ... % Gaussian
    + x(4)*(freq-x(3)) ... % linear baseline
    + x(5); % constant baseline


%%%%%%%%%%%%%%%% DOUBLE GAUSS MODEL %%%%%%%%%%%%%%%%
function F = DoubleGaussModel(x,freq)
% Function for DoubleGaussModel Model

% Two Gaussians
%  x(1) = gaussian amplitude 1
%  x(2) = width 1 ( 1/(2*sigma^2) )
%  x(3) = centre freq of peak 1
%  x(4) = gaussian amplitude 2
%  x(5) = width 2 ( 1/(2*sigma^2) )
%  x(6) = centre freq of peak 2
%  x(7) = amplitude of linear baseline
%  x(8) = constant amplitude offset

% MM: Allowing peaks to vary individually seems to work better than keeping
% the distance fixed (i.e., including J in the function)

F = x(1)*exp(x(2)*(freq-x(3)).*(freq-x(3))) + ...
    x(4)*exp(x(5)*(freq-x(6)).*(freq-x(6))) + ...
    x(7)*(freq-x(3)) + x(8);


%%%%%%%%%%%%%%%% TRIPLE GAUSS MODEL %%%%%%%%%%%%%%%%
function F = GABAGlxModel(x,freq)
% Function for GABA+Glx model

% Three Gaussians
%  x(1) = gaussian amplitude 1
%  x(2) = width 1 ( 1/(2*sigma^2) )
%  x(3) = centre freq of peak 1
%  x(4) = gaussian amplitude 2
%  x(5) = width 2 ( 1/(2*sigma^2) )
%  x(6) = centre freq of peak 2
%  x(7) = gaussian amplitude 3
%  x(8) = width 3 ( 1/(2*sigma^2) )
%  x(9) = centre freq of peak 3
%  x(10) = linear baseline slope
%  x(11) = sine baseline term
%  x(12) = cosine baseline term

% MM: Allowing peaks to vary individually seems to work better than keeping
% the distance fixed (i.e., including J in the function)

F = x(1)*exp(x(2)*(freq-x(3)).*(freq-x(3))) + ...
    x(4)*exp(x(5)*(freq-x(6)).*(freq-x(6))) + ...
    x(7)*exp(x(8)*(freq-x(9)).*(freq-x(9))) + ...
    x(10)*(freq-x(3)) + ...
    x(11)*sin(pi*freq/1.31/4) + x(12)*cos(pi*freq/1.31/4);


%%%%%%%%%%%%%%%% DOUBLE LORENTZ MODEL FOR ETOH %%%%%%%%%%%%%%%%
function F = EtOHModel(x,freq)
% Function for EtOH model

L1 = x(1) ./ (1 + ((freq-x(2)) / (x(3)/2)).^2);
L2 = x(4) ./ (1 + ((freq-x(5)) / (x(6)/2)).^2);
B  = x(7) .* (freq-x(3)) + x(8);
F  = L1 + L2 + B;


%%%%%%%%%%%%%%%% BASELINE %%%%%%%%%%%%%%%%
function F = BaselineModel(x,freq)
% Function for Baseline Model

F = x(2)*(freq-x(1))+x(3);


%%%%%%%%%%%%%%%% CALCULATE I.U. %%%%%%%%%%%%%%%%
function MRS_struct = CalcIU(MRS_struct, vox, metab, ii)
% Function for quantifying concentration in institutional units
% Convert metabolits and water areas to institutional units
% (pseudo-concentration in mmol/L)

TR = MRS_struct.p.TR(ii)/1e3;
TE = MRS_struct.p.TE(ii)/1e3;
if isfield(MRS_struct.p,'TR_water')
    TR_water = MRS_struct.p.TR_water(ii)/1e3;
else
    TR_water = TR;
end
if isfield(MRS_struct.p,'TE_water')
    TE_water = MRS_struct.p.TE_water(ii)/1e3;
else
    TE_water = TE;
end
PureWaterConc = 55000; % mmol/L
WaterVisibility = 0.65; % this is approx the value from Ernst, Kreis, Ross (1993, JMR)
T1_Water = 1.100; % average of WM and GM, Wansapura et al. 1999 (JMRI)
T2_Water = 0.095; % average of WM and GM, Wansapura et al. 1999 (JMRI)
N_H_Water = 2;

switch metab
    case 'GABA'
        EditingEfficiency = 0.5; % For TE = 68 ms
        T1_Metab = 1.31;  % Puts et al. 2013 (JMRI)
        T2_Metab = 0.088; % Edden et al. 2012 (JMRI)
        N_H_Metab = 2;
        MM = 0.45; % MM correction: fraction of GABA in GABA+ peak. (In TrypDep, 30 subjects: 55% of GABA+ was MM)
        % This fraction is platform and implementation dependent, based on length and
        % shape of editing pulses and ifis Henry method
        
    case 'Glx'
        EditingEfficiency = 0.4; % determined by FID-A simulations (for TE = 68 ms)
        T1_Metab = 1.23; % Posse et al. 2007 (MRM)
        T2_Metab = 0.18; % Ganji et al. 2012 (NMR Biomed)
        N_H_Metab = 1;
        MM = 1;
        
    case 'GSH'
        EditingEfficiency = 0.74;  % At 3T based on Quantification of Glutathione in the Human Brain by MR Spectroscopy at 3 Tesla:
        % Comparison of PRESS and MEGA-PRESS
        % Faezeh Sanaei Nezhad etal. DOI 10.1002/mrm.26532, 2016 -- MGSaleh
        T1_Metab = 0.40; % At 3T based on Doubly selective multiple quantum chemical shift imaging and
        % T1 relaxation time measurement of glutathione (GSH) in the human brain in vivo
        % In-Young Choi et al. NMR Biomed. 2013; 26: 28?34 -- MGSaleh
        T2_Metab = 0.12; % At 3T based on the ISMRM abstract
        % T2 relaxation times of 18 brain metabolites determined in 83 healthy volunteers in vivo
        % Milan Scheidegger et al. Proc. Intl. Soc. Mag. Reson. Med. 22 (2014)-- MGSaleh
        N_H_Metab = 2;
        MM = 1;
        
    case 'Lac'
        EditingEfficiency = 0.94; % determined by FID-A simulations (for TE = 140 ms)
        T1_Metab = 1.50; % Wijnen et al. 2015 (NMR Biomed)
        T2_Metab = 0.24; % Madan et al. 2015 (MRM) (NB: this was estimated in brain tumors)
        N_H_Metab = 3;
        MM = 1;
        
    case 'EtOH'
        EditingEfficiency = 0.5; % assuming same as GABA for now
        T1_Metab = 1.31;  % assuming same as GABA
        T2_Metab = 0.088; % assuming same as GABA
        N_H_Metab = 3;
        MM = 1;
end

T1_Factor = (1-exp(-TR_water./T1_Water)) ./ (1-exp(-TR./T1_Metab));
T2_Factor = exp(-TE_water./T2_Water) ./ exp(-TE./T2_Metab);

if strcmpi(MRS_struct.p.vendor,'Siemens_rda')
    % Factor of 2 is appropriate for averaged Siemens data (read in separately as ON and OFF)
    MRS_struct.out.(vox).(metab).ConcIU(ii) = (MRS_struct.out.(vox).(metab).Area(ii) ./ MRS_struct.out.(vox).water.Area(ii))  ...
        .* PureWaterConc .* WaterVisibility .* T1_Factor .* T2_Factor .* (N_H_Water./N_H_Metab) ...
        .* MM ./ 2 ./ EditingEfficiency;
else
    MRS_struct.out.(vox).(metab).ConcIU(ii) = (MRS_struct.out.(vox).(metab).Area(ii) ./ MRS_struct.out.(vox).water.Area(ii))  ...
        .* PureWaterConc .* WaterVisibility .* T1_Factor .* T2_Factor .* (N_H_Water./N_H_Metab) ...
        .* MM ./ EditingEfficiency;
end


%%%%%%%%%%%%%%%% INSET FIGURE %%%%%%%%%%%%%%%%
function [h_main, h_inset] = inset(main_handle, inset_handle,inset_size)
% Function for figure settings

% The function plotting figure inside figure (main and inset) from 2 existing figures.
% inset_size is the fraction of inset-figure size, default value is 0.35
% The outputs are the axes-handles of both.
%
% An examle can found in the file: inset_example.m
%
% Moshe Lindner, August 2010 (C).

if nargin == 2
    inset_size = 0.35;
end

inset_size = inset_size*.5;
new_fig = gcf;
main_fig = findobj(main_handle,'Type','axes');
h_main = copyobj(main_fig,new_fig);
set(h_main,'Position',get(main_fig,'Position'))
inset_fig = findobj(inset_handle,'Type','axes');
h_inset = copyobj(inset_fig,new_fig);
ax = get(main_fig,'Position');
set(h_inset,'Position', [1.3*ax(1)+ax(3)-inset_size, 1.001*ax(2)+ax(4)-inset_size, inset_size*0.7, inset_size*0.9])



