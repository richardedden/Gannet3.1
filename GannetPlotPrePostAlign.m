function GannetPlotPrePostAlign(MRS_struct, vox, ii, kk)
% Plot pre-/post-alignment spectra
% Updates by MGSaleh 2016, MM 2017-2019

if MRS_struct.p.HERMES
    
    SpectraToPlot = zeros(2*length(MRS_struct.p.target), length(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).diff));
    count = 0;
    for jj = 1:length(MRS_struct.p.target)
        SpectraToPlot((1:2)+count,:) = [MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:); ...
                                        MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:)];
        count = count + 2;
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
                if MRS_struct.p.phantom
                    peakrange(jj,:) = MRS_struct.spec.freq <= 4.25 & MRS_struct.spec.freq >= 1.0;
                else
                    peakrange(jj,:) = MRS_struct.spec.freq <= 4.1 & MRS_struct.spec.freq >= 2.26;
                end
            case 'GSH'
                if ~MRS_struct.p.HERCULES
                    peakrange(jj,:) = MRS_struct.spec.freq <= 3.5 & MRS_struct.spec.freq >= 0.5;
                else
                    peakrange(jj,:) = MRS_struct.spec.freq <= 4.25 & MRS_struct.spec.freq >= 3;
                end
            case {'Lac','EtOH'}
                peakrange(jj,:) = MRS_struct.spec.freq <= 3 & MRS_struct.spec.freq >= 0.5;
        end
    end
    
    specbaseline = mean(real(SpectraToPlot(:,baserange)),2);
    SpectraToPlot = SpectraToPlot - repmat(specbaseline, [1 size(SpectraToPlot,2)]);
    
    % Stack spectra
    if MRS_struct.p.phantom
        signalrange = max(max(real(SpectraToPlot(1:2,logical(peakrange(1,:)))))) - min(min(real(SpectraToPlot(1:2,logical(peakrange(1,:))))));
        SpectraToPlot(3:4,:) = SpectraToPlot(3:4,:) + signalrange;
        
        yaxismax = max(real(SpectraToPlot(3,logical(peakrange(2,:)))));
        yaxismin = min(real(SpectraToPlot(1,logical(peakrange(1,:)))));
        yrange = abs(yaxismax - yaxismin);
        yaxismax = yaxismax + 0.1*yrange;
        yaxismin = yaxismin - 0.1*yrange;
    else
        if all(ismember(MRS_struct.p.target,{'GABAGlx','GSH'})) || all(ismember(MRS_struct.p.target,{'GABA','GSH'}))
            if ~MRS_struct.p.HERCULES
                peakheight = abs(min(real(SpectraToPlot(3,logical(peakrange(2,:))))));
            else
                peakheight = abs(max(real(SpectraToPlot(1,logical(peakrange(1,:))))));
            end
            SpectraToPlot(3:4,:) = SpectraToPlot(3:4,:) + 1.5*peakheight;
            
            yaxismax = max(real(SpectraToPlot(3,logical(peakrange(2,:)))));
            yaxismin = min(real(SpectraToPlot(1,logical(peakrange(1,:)))));
            yrange = abs(yaxismax - yaxismin);
            yaxismax = yaxismax + 0.1*yrange;
            yaxismin = yaxismin - 0.1*yrange;
        elseif all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
            peakheight = abs(min(real(SpectraToPlot(5,logical(peakrange(3,:))))));
            SpectraToPlot(5:6,:) = SpectraToPlot(5:6,:) + 1.5*peakheight;
            peakheight = abs(max(real(SpectraToPlot(1,logical(peakrange(1,:))))));
            SpectraToPlot(1:2,:) = SpectraToPlot(1:2,:) - 1.5*peakheight;
            
            yaxismax = abs(max(real(SpectraToPlot(5,logical(peakrange(3,:))))));
            yaxismax = yaxismax + 0.2*yaxismax;
            yaxismin = max(real(SpectraToPlot(1,logical(peakrange(1,:)))));
            yaxismin = yaxismin - 4*abs(yaxismin);
        else
            peakheight = abs(min(real(SpectraToPlot(5,logical(peakrange(3,:))))));
            SpectraToPlot(5:6,:) = SpectraToPlot(5:6,:) + 1.5*peakheight;
        end
    end
    
    count = 0;
    for jj = 1:length(MRS_struct.p.target)
        hold on;
        plot(MRS_struct.spec.freq, real(SpectraToPlot(2+count*2,:)), 'Color', 'r');
        plot(MRS_struct.spec.freq, real(SpectraToPlot(1+count*2,:)), 'Color', 'b');
        hold off;
        count = count + 1;
    end
    
else
    
    SpectraToPlot = [MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).diff(ii,:); ...
                     MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).diff_noalign(ii,:)];
    
    % Shift baselines to zero
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
    
    specbaseline = mean(real(SpectraToPlot(:,baserange)),2);
    SpectraToPlot = SpectraToPlot - repmat(specbaseline, [1 length(SpectraToPlot)]);
    
    % Stack spectra
    peakheight = max(abs(real(SpectraToPlot(1,peakrange))));
    SpectraToPlot(2,:) = SpectraToPlot(2,:) + peakheight + 0.05*peakheight;
    hold on;
    plot(MRS_struct.spec.freq, real(SpectraToPlot(2,:)), 'Color', 'r');
    plot(MRS_struct.spec.freq, real(SpectraToPlot(1,:)), 'Color', 'b');
    hold off;
    
    yaxismax = max(max(real(SpectraToPlot(:,peakrange)),[],2));
    yaxismin = min(min(real(SpectraToPlot(:,peakrange)),[],2));
    yrange = abs(yaxismax - yaxismin);
    yaxismax = yaxismax + 0.1*yrange;
    if any(strcmp(MRS_struct.p.target{1}, {'GABA','Glx','GABAGlx'}))
        yaxismin = yaxismin - 0.3*yrange;
    else
        yaxismin = yaxismin - 0.1*yrange;
    end
    
end

if strcmp(MRS_struct.p.target{1},'EtOH')
    legend({'pre','post'},'EdgeColor',[1 1 1],'Location','northwest');
else
    legend({'pre','post'},'EdgeColor',[1 1 1]);
end
axis([0 5 yaxismin yaxismax]);
set(gca,'XDir','reverse','TickDir','out','box','off','XTick',0:5);
ax = get(gca,'YAxis');
set(ax,'Visible','off');
xlabel('ppm');
if MRS_struct.p.HERMES
    title({'Edited Spectra';'(pre- and post-alignment)'});
else
    title({'Edited Spectrum';'(pre- and post-alignment)'});
end
        


