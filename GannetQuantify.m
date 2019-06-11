function MRS_struct = GannetQuantify(MRS_struct)

MRS_struct.version.quantify = '190607';

% ******
% RE (190107): Major change to water concentration calc to bring into line
% with Gasparovic et al. 2006. Tagged '% Gasparovic et al. method (RE)'
% Fractions implemented as molar fractions, not volume fractions
% ******

cWM = 1; % concentration of GABA in pure WM
cGM = 2; % concentration of GABA in pure GM
alpha = cWM/cGM;

% Constants
% From Wansapura et al. 1999 (JMRI)
%        T1          T2
% WM   832 +/- 10  79.2 +/- 0.6
% GM  1331 +/- 13  110 +/- 2
%
% From Lu et al. 2005 (JMRI)
% CSF T1 = 3817 +/- 424msec - but state may underestimated and that 4300ms
% is likely more accurate - but the reference is to an ISMRM 2001 abstract
% MacKay (last author) 2006 ISMRM abstract has T1 CSF = 3300 ms
% CSF T2 = 503.0 +/- 64.3 Piechnik MRM 2009; 61: 579
% However, other values from Stanisz et al:
% CPMG for T2, IR for T1
% T2GM = 99 +/ 7, lit: 71+/- 10 (27)
% T1GM = 1820 +/- 114, lit 1470 +/- 50 (29)
% T2WM = 69 +/-3 lit 56 +/- 4 (27)
% T1WM = 1084 +/- 45 lit 1110 +/- 45 (29)

T1w_WM    = 0.832;
T2w_WM    = 0.0792;
T1w_GM    = 1.331;
T2w_GM    = 0.110;
T1w_CSF   = 3.817;
T2w_CSF   = 0.503;
N_H_Water = 2;

% Determine concentration of water in GM, WM and CSF
% Gasparovic et al. 2006 (MRM) uses relative densities, ref to
% Ernst et al. 1993 (JMR)
% fGM = 0.78
% fWM = 0.65
% fCSF = 0.97
% such that
% concw_GM = 0.78 * 55.51 mol/kg = 43.30
% concw_WM = 0.65 * 55.51 mol/kg = 36.08
% concw_CSF = 0.97 * 55.51 mol/kg = 53.84

concW_GM    = 43.30*1e3;
concW_WM    = 36.08*1e3;
concW_CSF   = 53.84*1e3;
molal_concW = 55.51*1e3; % Gasparovic et al. method (RE)

if MRS_struct.p.PRIAM
    vox = MRS_struct.p.Vox;
else
    vox = MRS_struct.p.Vox(1);
end

numscans = MRS_struct.p.numscans;

for ii = 1:numscans
    
    target = MRS_struct.p.target;
    
    tmp = strcmp(target,'GABAGlx');
    if any(tmp)
        if MRS_struct.p.HERMES
            target = {'GABA','Glx',target{~tmp}};
        else
            target = {'GABA','Glx'};
        end
    end
    
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
    
    % Loop over voxels if PRIAM
    for kk = 1:length(vox)
        
        meanfGM = mean(MRS_struct.out.(vox{kk}).tissue.fGM); % average GM fraction across subjects
        meanfWM = mean(MRS_struct.out.(vox{kk}).tissue.fWM); % average WM fraction across subjects
        
        fGM  = MRS_struct.out.(vox{kk}).tissue.fGM(ii);
        fWM  = MRS_struct.out.(vox{kk}).tissue.fWM(ii);
        fCSF = MRS_struct.out.(vox{kk}).tissue.fCSF(ii);
        
        % Gasparovic et al. method (RE)
        % Calculate molal fractions from volume fractions (equivalent to eqs. 5-7 in Gasparovic et al., 2006)
        molal_fGM  = (fGM*concW_GM) / (fGM*concW_GM + fWM*concW_WM + fCSF*concW_CSF);
        molal_fWM  = (fWM*concW_WM) / (fGM*concW_GM + fWM*concW_WM + fCSF*concW_CSF);
        molal_fCSF = (fCSF*concW_CSF) / (fGM*concW_GM + fWM*concW_WM + fCSF*concW_CSF);
        
        CorrFactor = (meanfGM + alpha*meanfWM) / ((fGM + alpha*fWM) * (meanfGM + meanfWM));
        
        for jj = 1:length(target)
            
            switch target{jj}
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
                    T1_Metab = 0.40 ; % At 3T based on Doubly selective multiple quantum chemical shift imaging and
                    % T1 relaxation time measurement of glutathione (GSH) in the human brain in vivo
                    % In-Young Choi et al. NMR Biomed. 2013; 26: 28-34 -- MGSaleh
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
            end
            
            % Gasparovic et al. method (RE)
            MRS_struct.out.(vox{kk}).(target{jj}).ConcIU_TissCorr(ii) = ...
                (MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) / MRS_struct.out.(vox{kk}).water.Area(ii)) * ...
                (N_H_Water / N_H_Metab) * MM / EditingEfficiency * molal_concW * ...
                (molal_fGM  * (1 - exp(-TR_water/T1w_GM)) * exp(-TE_water/T2w_GM) / ((1 - exp(-TR/T1_Metab)) * exp(-TE/T2_Metab)) + ...
                molal_fWM  * (1 - exp(-TR_water/T1w_WM)) * exp(-TE_water/T2w_WM) / ((1 - exp(-TR/T1_Metab)) * exp(-TE/T2_Metab)) + ...
                molal_fCSF * (1 - exp(-TR_water/T1w_CSF)) * exp(-TE_water/T2w_CSF) / ((1 - exp(-TR/T1_Metab)) * exp(-TE/T2_Metab))) / ...
                (1 - molal_fCSF);
            
            ConcIU_TissCorr_Harris = ...
                (MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) / MRS_struct.out.(vox{kk}).water.Area(ii)) * ...
                (N_H_Water / N_H_Metab) * MM / EditingEfficiency * ...
                (fGM * concW_GM * (1 - exp(-TR_water/T1w_GM)) * exp(-TE_water/T2w_GM) / ((1 - exp(-TR/T1_Metab)) * exp(-TE/T2_Metab)) + ...
                fWM * concW_WM * (1 - exp(-TR_water/T1w_WM)) * exp(-TE_water/T2w_WM) / ((1 - exp(-TR/T1_Metab)) * exp(-TE/T2_Metab)) + ...
                fCSF * concW_CSF * (1 - exp(-TR_water/T1w_CSF)) * exp(-TE_water/T2w_CSF) / ((1 - exp(-TR/T1_Metab)) * exp(-TE/T2_Metab)));
            MRS_struct.out.(vox{kk}).(target{jj}).ConcIU_AlphaTissCorr(ii) = ConcIU_TissCorr_Harris / (fGM + alpha*fWM);
            MRS_struct.out.(vox{kk}).(target{jj}).ConcIU_AlphaTissCorr_GrpNorm(ii) = ConcIU_TissCorr_Harris * CorrFactor;
            
        end
        
        % Build output figure
        if ishandle(105)
            clf(105);
        end
        h = figure(105);
        % Open figure in center of screen
        scr_sz = get(0, 'ScreenSize');
        fig_w = 1000;
        fig_h = 707;
        set(h,'Position',[(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
        set(h,'Color',[1 1 1]);
        figTitle = 'GannetQuantify Output';
        set(h,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');
        
        % Voxel co-registration        
        img_montage = MRS_struct.mask.(vox{kk}).img_montage{ii};
        img_montage = [img_montage(:,1:size(img_montage,2)/2); img_montage(:,size(img_montage,2)/2+1:end)];
        
        ha = subplot(2,2,1);
        imagesc(img_montage);
        colormap('gray');
        img = MRS_struct.mask.(vox{kk}).img{ii}(:);
        caxis([0 mean(img(img > 0.01)) + 3*std(img(img > 0.01))]);
        axis equal;
        axis tight;
        axis off;
        pos = get(ha,'pos');
        s = 0.04;
        set(ha,'pos',[pos(1)-s, pos(2)-s-0.02, pos(3)+2*s, pos(4)+2*s]);
        
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
            [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii*2-1});
        else
            [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii});
        end
        fname = [tmp tmp2];
        if length(fname) > 30
            fname = [fname(1:30) '...'];
        end
        [~,tmp3,tmp4] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
        T1image = [tmp3 tmp4];
        if length(T1image) > 30
            T1image = [T1image(1:30) '...'];
        end
        title(['Voxel from ' fname ' on ' T1image], 'Interpreter', 'none');
        
        % Post-alignment spectra + model fits
        subplot(2,2,3);
        GannetPlotPrePostAlign2(MRS_struct, vox, ii);
        
        % Output results
        subplot(2,2,2);
        axis off;
        
        target = MRS_struct.p.target;
        
        tmp = strcmp(target,'GABAGlx');
        if any(tmp)
            if MRS_struct.p.HERMES
                target = {'GABA','Glx',target{~tmp}};
            else
                target = {'GABA','Glx'};
            end
        end
        
        for jj = 1:length(target)
            
            switch target{jj}
                case 'GABA'
                    tmp2 = 'GABA+';
                case {'Glx','GSH','Lac'}
                    tmp2 = target{jj};
            end
            
            shift = 0;
            for ll = 1:3
                text_pos = 0.9;
                if ll == 1
                    tmp1 = 'Relaxation-, tissue-corrected (Gasparovic et al. method)';
                    tmp3 = sprintf('%.2f i.u.', MRS_struct.out.(vox{kk}).(target{jj}).ConcIU_TissCorr(ii));
                    text(0, text_pos, tmp1, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 10);
                elseif ll == 2
                    text_pos = text_pos - 0.2 - shift;
                    tmp1 = 'Relaxation-, tissue-, alpha-corrected (Harris et al. method)';
                    tmp3 = sprintf('%.2f i.u.', MRS_struct.out.(vox{kk}).(target{jj}).ConcIU_AlphaTissCorr(ii));
                    text(0, text_pos, tmp1, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 10);
                elseif ll == 3
                    text_pos = text_pos - 0.4 - shift;
                    tmp1a = 'Relaxation-, tissue-, alpha-corrected; average-voxel-normalized';
                    tmp1b = '(Harris et al. method)';
                    tmp3 = sprintf('%.2f i.u.', MRS_struct.out.(vox{kk}).(target{jj}).ConcIU_AlphaTissCorr_GrpNorm(ii));
                    text(0, text_pos, tmp1a, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 10);
                    text(0, text_pos - 0.1, tmp1b, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 10);
                end
                if ll == 3
                    text_pos = text_pos - 0.1*(jj+1);
                else
                    text_pos = text_pos - 0.1*jj;
                end
                text(0.4, text_pos, [tmp2 '/Water: '], 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                text(0.425, text_pos, tmp3, 'FontName', 'Arial', 'FontSize', 10);
                if MRS_struct.p.HERMES
                    shift = shift + 0.1*(numel(target)-1);
                else
                    shift = shift + 0.1;
                end
            end
            
        end
        
        tmp1 = 'Filename:';
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
            [~,tmp2,tmp3] = fileparts(MRS_struct.metabfile{ii*2-1});
        else
            [~,tmp2,tmp3] = fileparts(MRS_struct.metabfile{ii});
        end
        fname = [tmp2 tmp3];
        if length(fname) > 30
            fname = [fname(1:30) '...'];
        end
        text(0.4, text_pos-0.1, tmp1, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
        text(0.425, text_pos-0.1, fname, 'FontName', 'Arial', 'FontSize', 10, 'Interpreter', 'none');
        
        tmp1 = 'Anatomical image:';
        [~,tmp2,tmp3] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
        T1image = [tmp2 tmp3];
        if length(T1image) > 30
            T1image = [T1image(1:30) '...'];
        end
        text(0.4, text_pos-0.2, tmp1, 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
        text(0.425, text_pos-0.2, T1image, 'FontName', 'Arial', 'FontSize', 10, 'Interpreter', 'none');
        
        text(0.4, text_pos-0.3, 'QuantifyVer:', 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
        text(0.425, text_pos-0.3, MRS_struct.version.quantify, 'FontName', 'Arial', 'FontSize', 10);
        
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
        if ~exist(fullfile(pwd, 'GannetQuantify_output'),'dir')
            mkdir(fullfile(pwd, 'GannetQuantify_output'));
        end
        
        % Save PDF output
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[11 8.5]);
        set(h,'PaperPosition',[0 0 11 8.5]);
        
        if strcmpi(MRS_struct.p.vendor,'Philips_data')
            pdfname = fullfile(pwd, 'GannetQuantify_output', [fullpath '_' vox{kk} '_quantify.pdf']);
        else
            pdfname = fullfile(pwd, 'GannetQuantify_output', [metabfile_nopath '_' vox{kk} '_quantify.pdf']);
        end
        saveas(h, pdfname);
        
        if ii == numscans
            if MRS_struct.p.mat % save MRS_struct as mat file
                mat_name = ['GannetQuantify_output/MRS_struct_' vox{kk} '.mat'];
                save(mat_name,'MRS_struct');
            end
            if MRS_struct.p.csv % export MRS_struct fields into csv file
                ExportToCSV(MRS_struct, MRS_struct.p.target, kk, 'quantify');
            end
        end
        
    end
    
end

end



