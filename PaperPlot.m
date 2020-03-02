function PaperPlot(MRS_struct, varargin)
% PaperPlot(MRS_struct, varargin)
%
% This function will plot the difference spectra saved in MRS_struct. The
% corresponding model fits optionally can also be plotted. Users can choose
% to plot a single spectrum, a select number of spectra or all spectra.
% Multiple spectra will be overlaid in the same figure. If data were
% acquired with HERMES, then each Hadamard-combined difference spectrum
% will be plotted in separate subplots.
%
% To export plots at publication quality, consider using PaperPlot with
% Yair Altman's excellent export_fig toolbox
% (https://github.com/altmany/export_fig).
%
% Inputs:
%   MRS_struct: Structure output from GannetFit (required).
%   varargin: Optional inputs (entered as parameter-value pairs).
%           target:     (For HERMES data only.) Choose a single target
%                       metabolite to plot, entered as a string. Default is
%                       plotting of difference spectra for all target
%                       metabolites.
%           specNum:    Spectra to plot, entered as a scalar or vector. All
%                       spectra are plotted by default.
%           freqLim:    Limits of ppm axis, entered as a two-element
%                       vector. Default is [0.5 4.5].
%           signalLim:  Limits of signal axis, entered as a two-element
%                       vector. Default is an empty vector (automatic
%                       scaling; recommended).
%           plotModel:  Plot signal model fit(s), entered as a logical.
%                       Default is false.
%           plotAvg:    Plot the group-average spectrum, entered as a
%                       logical. Default is false.
%           plotStd:    If plotAvg is true, also show the +/- 1 standard
%                       deviation, entered as a logical. Default is false.
%           plotCI:     If plotAvg is true, also show the 95% confidence
%                       interval, entered as a logical. Default is false.
%
% Examples:
%   PaperPlot(MRS_struct, 'specNum', [1 3 4]);
%       This will plot the 1st, 3rd and 4th difference spectra in
%       MRS_struct along with the model fits of the peak(s) specified in
%       MRS_struct.p.target.
%
%   PaperPlot(MRS_struct, 'freqLim', [2.5 3.5], 'plotModel', true);
%       This will plot all difference spectra in MRS_struct with the model
%       fits of the peak(s) and limit the ppm axis from 2.5 to 3.5 ppm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Parse inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    error('Not enough inputs! MRS_struct is required!');
end

% Set some defaults
vox = MRS_struct.p.Vox;
if ~MRS_struct.p.PRIAM
    vox = vox(1);
end
defaultTarget = MRS_struct.p.target;
defaultnSpec = 1:length(MRS_struct.metabfile);
defaultFreqLim = [0.5 4.5];
defaultSignalLim = [];
defaultPlotModel = false;
defaultPlotAvg = false;
defaultPlotStd = false;
defaultPlotCI = false;
expectedTargets = {'GABAGlx','GSH','Lac','EtOH','GABA','Glx'};
grey = [0.6 0.6 0.6];
shading = 0.3;

% Parse input arguments
p = inputParser;
p.CaseSensitive = false;
p.addParameter('target', defaultTarget, @(x) any(validatestring(x,expectedTargets)));
p.addParameter('specNum', defaultnSpec);
p.addParameter('freqLim', defaultFreqLim);
p.addParameter('signalLim', defaultSignalLim);
p.addParameter('plotModel', defaultPlotModel, @(x) islogical(x));
p.addParameter('plotAvg', defaultPlotAvg, @(x) islogical(x));
p.addParameter('plotStd', defaultPlotStd, @(x) islogical(x));
p.addParameter('plotCI', defaultPlotCI, @(x) islogical(x));
p.parse(varargin{:});

target = p.Results.target;
if ischar(target)
    target = {target};
end
specNum = p.Results.specNum;
freqLim = p.Results.freqLim;
signalLim = p.Results.signalLim;
plotModel = p.Results.plotModel;
plotAvg = p.Results.plotAvg;
plotStd = p.Results.plotStd;
plotCI = p.Results.plotCI;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Plot spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq = MRS_struct.spec.freq;

if isfield(MRS_struct.spec.(vox{1}).(target{1}),'diff_scaled')
    diff = 'diff_scaled';
    fprintf('\nNB: Spectra are normalized to the amplitude of the respective modeled unsuppressed water reference signal.\n\n');
else
    diff = 'diff';
end

for ii = 1:length(vox)
    
    H = figure(199+ii);
    scr_sz = get(0, 'ScreenSize');
    fig_w = 1000;
    if length(target) > 1
        fig_h = 1000;
    else
        fig_h = 500;
    end
    set(H, 'Color', 'w', 'Position', [(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
    clf;
    
    if isfield(MRS_struct.out.(vox{ii}),'water')
        scaleFactor = MRS_struct.out.(vox{ii}).water.ModelParam(specNum,1);
    else
        scaleFactor = ones(1,length(specNum));
    end
    
    for jj = 1:length(target)
        
        if length(target) > 1
            H = subplot(length(target),1,jj);
        else
            H = gca;
        end
        
        switch target{jj}
            case 'GABA'
                if length(target) == 3 && all(ismember(target,{'EtOH','GABA','GSH'}))
                    modelFreq = freq(freq <= 3.55 & freq >= 2.6);
                else
                    modelFreq = freq(freq <= 3.55 & freq >= 2.79);
                end
                model = @GaussModel;
                baselineFreq = freq <= 3.5 & freq >= 3.4;
            case 'GABAGlx'
                modelFreq = freq(freq <= 4.1 & freq >= 2.79);
                model = @GABAGlxModel;
                baselineFreq = freq <= 3.5 & freq >= 3.4;
            case 'GSH'
                modelFreq = freq(freq <= 3.5 & freq >= 2.25);
                if MRS_struct.p.TE(1) < 100
                    model = @FiveGaussModel;
                else
                    model = @SixGaussModel;
                end
                baselineFreq = freq <= 1.8 & freq >= 1.7;
            case 'Lac'
                modelFreq = freq(freq <= 1.8 & freq >= 0.5);
                model = @FourGaussModel;
            case 'EtOH'
                modelFreq = freq(freq <= 1.8 & freq >= 0.6);
                model = @EtOHModel;
        end
        
        % Demean baseline
        baseMean = repmat(mean(real(MRS_struct.spec.(vox{ii}).(target{jj}).(diff)(specNum,baselineFreq)),2), ...
            [1 size(MRS_struct.spec.(vox{ii}).(target{jj}).(diff),2)]);
        
        if numel(specNum) > 1 && plotAvg
            
            % Find mean, std and 95% CI
            mu = mean(real(MRS_struct.spec.(vox{ii}).(target{jj}).(diff)(specNum,:)) - baseMean,1);
            sigma = std(real(MRS_struct.spec.(vox{ii}).(target{jj}).(diff)(specNum,:)) - baseMean,[],1);
            stderr = sigma / sqrt(numel(specNum));
            UB.sigma = mu + sigma;
            LB.sigma = mu - sigma;
            UB.ci = mu + 1.96 * stderr;
            LB.ci = mu - 1.96 * stderr;
            
            hold on;
            if plotStd && plotCI
                patch([freq fliplr(freq)], [real(UB.sigma) fliplr(real(LB.sigma))], 1, 'FaceColor', grey+(1-grey)*(1-shading), 'EdgeColor', 'none');
                patch([freq fliplr(freq)], [real(UB.ci) fliplr(real(LB.ci))], 1, 'FaceColor', (grey-0.4)+(1-(grey-0.4))*(1-shading), 'EdgeColor', 'none');
            elseif plotStd
                patch([freq fliplr(freq)], [real(UB.sigma) fliplr(real(LB.sigma))], 1, 'FaceColor', grey+(1-grey)*(1-shading), 'EdgeColor', 'none');
            elseif plotCI
                patch([freq fliplr(freq)], [real(UB.ci) fliplr(real(LB.ci))], 1, 'FaceColor', grey+(1-grey)*(1-shading), 'EdgeColor', 'none');
            end
            h = plot(freq, mu, 'k', 'LineWidth', 1);
            hold off;
            
        else
            
            hold on;
            for kk = 1:numel(specNum)
                if plotModel
                    if strcmp(target{jj},'GABAGlx')
                        h(:,kk) = plot(freq, real(MRS_struct.spec.(vox{ii}).(target{jj}).(diff)(specNum(kk),:)) - baseMean(kk,:), 'k', ...
                            modelFreq, model(MRS_struct.out.(vox{ii}).GABA.ModelParam(specNum(kk),:),modelFreq) ./ scaleFactor(kk) - baseMean(kk,1), 'r', 'LineWidth', 1);
                    else
                        h(:,kk) = plot(freq, real(MRS_struct.spec.(vox{ii}).(target{jj}).(diff)(specNum(kk),:)) - baseMean(kk,:), 'k', ...
                            modelFreq, model(MRS_struct.out.(vox{ii}).(target{jj}).ModelParam(specNum(kk),:),modelFreq) ./ scaleFactor(kk) - baseMean(kk,1), 'r', 'LineWidth', 1);
                    end
                else
                    h(:,kk) = plot(freq, real(MRS_struct.spec.(vox{ii}).(target{jj}).(diff)(specNum(kk),:)) - baseMean(kk,:), 'k', 'LineWidth', 1);
                end
            end
            hold off;
            
        end
        
        set(gca,'TickDir','out','XLim',freqLim,'XDir','reverse','Box','off','FontSize',20,'LineWidth',1);
        ax = get(gca,'YAxis');
        set(ax,'Visible','off');
        xlabel('ppm','FontWeight','bold','FontSize',28);
        
        % Set YLim
        if isempty(signalLim)
            switch target{jj}
                case {'GABAGlx','GABA','Glx'}
                    peakrange = freq <= 4.1 & freq >= 2.26;
                case 'GSH'
                    peakrange = freq <= 3.5 & freq >= 0.5;
                case {'Lac','EtOH'}
                    peakrange = freq <= 3 & freq >= 0.5;
            end
            if plotStd && plotCI
                yaxismax = max([UB.sigma(peakrange) UB.ci(peakrange)]);
                yaxismin = min([LB.sigma(peakrange) LB.ci(peakrange)]);
            elseif plotStd
                yaxismax = max(UB.sigma(peakrange));
                yaxismin = min(LB.sigma(peakrange));
            elseif plotCI
                yaxismax = max(UB.ci(peakrange));
                yaxismin = min(LB.ci(peakrange));
            else
                for kk = 1:size(h,2)
                    maxPeakHeight(kk) = max(h(1,kk).YData(peakrange)); %#ok<*AGROW>
                    minPeakHeight(kk) = min(h(1,kk).YData(peakrange));
                end
                yaxismax = max(maxPeakHeight);
                yaxismin = min(minPeakHeight);
            end
            yrange = abs(yaxismax - yaxismin);
            yaxismax = yaxismax + 0.1*yrange;
            if any(strcmp(target{jj}, {'GABAGlx','GABA','Glx'}))
                yaxismin = yaxismin - 0.3*yrange;
            else
                yaxismin = yaxismin - 0.1*yrange;
            end
            signalLim = [yaxismin yaxismax];
        end
        set(H,'YLim',signalLim);
        signalLim = [];
        
    end
    
    set(findall(H,'-property','FontName'),'FontName','Arial');
    
end



