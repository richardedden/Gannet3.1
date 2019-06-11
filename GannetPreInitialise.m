function MRS_struct = GannetPreInitialise(MRS_struct)

% Acquisition parameters
    MRS_struct.p.target = {'GABAGlx'}; % edited metabolite(s) of interest; only possible options are:
                                       % if MEGA-PRESS:
                                       %   {'GABAGlx'}, {'GSH'}, {'Lac'} or {'EtOH'}
                                       % if HERMES:
                                       %   {'GABAGlx','GSH'}, {'Lac','GSH'} or {'EtOH','GABA','GSH'}
                                       % if HERCULES:
                                       %   {'GABAGlx','GSH'}
    MRS_struct.p.ONOFForder = 'offfirst'; % order of editing pulses; options are 'onfirst' or 'offfirst'
    MRS_struct.p.Water_Positive = 1; % 1 = yes; if water suppression method inverts residual water signal, set to 0
    MRS_struct.p.seqorig = 'JHU'; % origin of Philips patch; options are 'JHU' or 'Philips'
    
% Analysis parameters
    MRS_struct.p.LB = 3; % line-broadening (in Hz)
    MRS_struct.p.water_phase_correction = 1; % 1 = YES; perform phase correction
    MRS_struct.p.data_phase_correction = 0; % 1 = YES; perform phase correction
    MRS_struct.p.water_removal = 1; % 1 = YES; remove residual water using HSVD in HERMES (recommended for GABA/GSH editing) or GSH-edited MEGA-PRESS data
    MRS_struct.p.AlignTo = 'RobustSpecReg'; % options are 'RobustSpecReg' (recommended), 'SpecReg', 'SpecRegHERMES',
                                            % 'Cr', 'Cho', 'NAA', 'H2O', 'CrOFF' or 'none'
    MRS_struct.p.Vox = {'vox1','vox2'}; % for naming voxels acquired by PRIAM, e.g.: {'anterior','posterior'}, {'right','left'}, etc.
    MRS_struct.p.GSH_model = 'FiveGauss'; % choice of model for fitting GSH;
                                          % options are 'FiveGauss' (recommended for medium-TE HERMES) or 'SixGauss' (recommended for long-TE MEGA-PRESS)
    MRS_struct.p.FitResidWater = 0; % 1 = YES, fit the residual water peak in the DIFF spectrum to calculate water suppression factor
    
% Flags
    MRS_struct.p.HERMES   = 0; % 1 = YES, 0 = NO (for MEGA-PRESS)
    MRS_struct.p.HERCULES = 0; % 1 = YES, 0 = NO (for MEGA-PRESS; MRS_struct.p.HERMES should be set to 1 as well)
    MRS_struct.p.PRIAM    = 0; % 1 = YES, 0 = NO
    MRS_struct.p.phantom  = 0; % 1 = YES (assumes phantom was scanned at room temperature), 0 = NO (for in vivo data)
    MRS_struct.p.mat      = 0; % 1 = YES, save MRS_struct as .mat file
    MRS_struct.p.sdat     = 0; % 1 = YES, save processed difference spectrum as .sdat file (only for Philips SDAT MEGA-PRESS datasets)
    MRS_struct.p.csv      = 0; % 1 = YES, extract useful data from MRS_struct and export to .csv file (applies to GannetFit, GannetSegment and GannetQuantify)
    
end



