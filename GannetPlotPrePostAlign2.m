function GannetPlotPrePostAlign2(MRS_struct, vox, ii)
% Plot pre-/post-alignment spectra
% Updates by MGSaleh 2016, MM 2017-2019

for kk = 1:length(vox)
    
    if MRS_struct.p.HERMES
        
        SpectraToPlot = zeros(length(MRS_struct.p.target), length(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).diff));
        for jj = 1:length(MRS_struct.p.target)
            SpectraToPlot(jj,:) = MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:); %#ok<*AGROW>
        end
        
        model = cell(1,2);
        freqbounds = cell(1,2);
        for jj = 1:length(MRS_struct.p.target)
            switch MRS_struct.p.target{jj}
                case 'GABA'
                    freqbounds{jj} = MRS_struct.spec.freq <= 3.55 & MRS_struct.spec.freq >= 2.79;
                    model{jj} = GaussModel(MRS_struct.out.(vox{kk}).GABA.ModelParam(ii,:),MRS_struct.spec.freq(freqbounds{jj}));
                case 'GSH'
                    freqbounds{jj} = MRS_struct.spec.freq <= 3.3 & MRS_struct.spec.freq >= 2.35;
                    model{jj} = FiveGaussModel(MRS_struct.out.(vox{kk}).GSH.ModelParam(ii,:),MRS_struct.spec.freq(freqbounds{jj}));
                case 'GABAGlx'
                    freqbounds{jj} = MRS_struct.spec.freq <= 4.1 & MRS_struct.spec.freq >= 2.79;
                    model{jj} = GABAGlxModel(MRS_struct.out.(vox{kk}).GABA.ModelParam(ii,:),MRS_struct.spec.freq(freqbounds{jj}));
                case 'Lac'
                    freqbounds{jj} = MRS_struct.spec.freq <= 1.8 & MRS_struct.spec.freq >= 0.5;
                    model{jj} = FourGaussModel(MRS_struct.out.(vox{kk}).Lac.ModelParam(ii,:),MRS_struct.spec.freq(freqbounds{jj}));
            end
        end
        
        % Shift baselines to zero
        baserange = MRS_struct.spec.freq <= 0 & MRS_struct.spec.freq >= -0.5;
        % Some bandwidth-limited acquisitions may not record anything below
        % 0 ppm, in this case get the baseline from the other side of
        % water. (GO: 180213)
        if sum(baserange) == 0
            baserange = MRS_struct.spec.freq >= 7 & MRS_struct.spec.freq <= 8;
        end
        peakrange = zeros(length(MRS_struct.p.target), length(MRS_struct.spec.freq));
        for jj = 1:length(MRS_struct.p.target)
            switch MRS_struct.p.target{jj}
                case {'GABA','Glx','GABAGlx'}
                    peakrange(jj,:) = MRS_struct.spec.freq <= 4.1 & MRS_struct.spec.freq >= 2.26;
                case 'GSH'
                    peakrange(jj,:) = MRS_struct.spec.freq <= 3.5 & MRS_struct.spec.freq >= 0.5;
                case {'Lac','EtOH'}
                    peakrange(jj,:) = MRS_struct.spec.freq <= 3 & MRS_struct.spec.freq >= 0.5;
            end
        end
        
        specbaseline = mean(real(SpectraToPlot(:,baserange)),2);
        SpectraToPlot = SpectraToPlot - repmat(specbaseline, [1 length(SpectraToPlot)]);
        for jj = 1:length(MRS_struct.p.target)
            model{jj} = model{jj} - specbaseline(jj);
        end
        
        % Stack spectra (MM: 180724)
        if all(ismember(MRS_struct.p.target,{'GABAGlx','GSH'})) || all(ismember(MRS_struct.p.target,{'GABA','GSH'}))
            peakheight = abs(min(real(SpectraToPlot(2,logical(peakrange(2,:))))));
            SpectraToPlot(2,:) = SpectraToPlot(2,:) + 1.5*peakheight;
            model{2} = model{2} + 1.5*peakheight;
            
            yaxismax = max(real(SpectraToPlot(2,logical(peakrange(2,:)))));
            yaxismin = min(real(SpectraToPlot(1,logical(peakrange(1,:)))));
            yrange = abs(yaxismax - yaxismin);
            yaxismax = yaxismax + 0.1*yrange;
            yaxismin = yaxismin - 0.1*yrange;
        elseif all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
            peakheight = abs(min(real(SpectraToPlot(3,logical(peakrange(3,:))))));
            SpectraToPlot(3,:) = SpectraToPlot(3,:) + 1.5*peakheight;
            peakheight = abs(max(real(SpectraToPlot(1,logical(peakrange(1,:))))));
            SpectraToPlot(1,:) = SpectraToPlot(1,:) - 1.5*peakheight;
            
            yaxismax = abs(max(real(SpectraToPlot(3,logical(peakrange(3,:))))));
            yaxismax = yaxismax + 0.2*yaxismax;
            yaxismin = max(real(SpectraToPlot(1,logical(peakrange(1,:)))));
            yaxismin = yaxismin - 4*abs(yaxismin);
        end
        
        for jj = 1:length(MRS_struct.p.target)
            hold on;
            plot(MRS_struct.spec.freq, real(SpectraToPlot(jj,:)), 'Color', 'k');
            plot(MRS_struct.spec.freq(freqbounds{jj}), model{jj}, 'Color', 'r');
            hold off;
        end
        
    else
        
        SpectraToPlot = MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).diff(ii,:);
        
        switch MRS_struct.p.target{1}
            case 'GABA'
                freqbounds = MRS_struct.spec.freq <= 3.55 & MRS_struct.spec.freq >= 2.79;
                model = GaussModel(MRS_struct.out.(vox{kk}).GABA.ModelParam(ii,:),MRS_struct.spec.freq(freqbounds));
            case 'GSH'
                freqbounds = MRS_struct.spec.freq <= 3.3 & MRS_struct.spec.freq >= 2.35;
                model = FiveGaussModel(MRS_struct.out.(vox{kk}).GSH.ModelParam(ii,:),MRS_struct.spec.freq(freqbounds));
            case 'GABAGlx'
                freqbounds = MRS_struct.spec.freq <= 4.1 & MRS_struct.spec.freq >= 2.79;
                model = GABAGlxModel(MRS_struct.out.(vox{kk}).GABA.ModelParam(ii,:),MRS_struct.spec.freq(freqbounds));
            case 'Lac'
                freqbounds = MRS_struct.spec.freq <= 1.8 & MRS_struct.spec.freq >= 0.5;
                model = FourGaussModel(MRS_struct.out.(vox{kk}).Lac.ModelParam(ii,:),MRS_struct.spec.freq(freqbounds));
        end
        
        % Shift baselines to zero (MM: 180108)
        baserange = MRS_struct.spec.freq <= 0 & MRS_struct.spec.freq >= -0.5;
        % Some bandwidth-limited acquisitions may not record anything below
        % 0 ppm, in this case get the baseline from the other side of
        % water. (GO: 180213)
        if sum(baserange) == 0
            baserange = MRS_struct.spec.freq >= 7 & MRS_struct.spec.freq <= 8;
        end
        switch MRS_struct.p.target{1}
            case {'GABA','Glx','GABAGlx'}
                peakrange = MRS_struct.spec.freq <= 4.1 & MRS_struct.spec.freq >= 2.26;
            case 'GSH'
                peakrange = MRS_struct.spec.freq <= 3.5 & MRS_struct.spec.freq >= 0.5;
            case {'Lac','EtOH'}
                peakrange = MRS_struct.spec.freq <= 3 & MRS_struct.spec.freq >= 0.5;
        end
        
        specbaseline = mean(real(SpectraToPlot(baserange)));
        SpectraToPlot = SpectraToPlot - specbaseline;
        model = model - specbaseline;
        
        hold on;
        plot(MRS_struct.spec.freq, real(SpectraToPlot), 'Color', 'k');
        plot(MRS_struct.spec.freq(freqbounds), model, 'Color', 'r');
        hold on;
        
        yaxismax = max(real(SpectraToPlot(peakrange)));
        yaxismin = min(real(SpectraToPlot(peakrange)));
        yrange = abs(yaxismax - yaxismin);
        yaxismax = yaxismax + 0.1*yrange;
        if any(strcmp(MRS_struct.p.target{1}, {'GABA','Glx','GABAGlx'}))
            yaxismin = yaxismin - 0.3*yrange;
        else
            yaxismin = yaxismin - 0.1*yrange;
        end
        
    end
    
    axis([0 5 yaxismin yaxismax]);
    set(gca,'XDir','reverse','TickDir','out','box','off','XTick',0:5);
    ax = get(gca,'YAxis');
    set(ax,'Visible','off');
    xlabel('ppm');
    if MRS_struct.p.HERMES
        title('Edited Spectra and Model Fits');
    else
        title('Edited Spectrum and Model Fit');
    end
    
end



