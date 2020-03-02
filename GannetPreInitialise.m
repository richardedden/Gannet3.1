function MRS_struct = GannetPreInitialise(MRS_struct)

% Acquisition parameters
    MRS_struct.p.target = {'GABAGlx'}; % edited metabolite(s) of interest; allowable options are:
                                       % if MEGA-PRESS:
                                       %   {'GABAGlx'}, {'GSH'}, {'Lac'}, or {'EtOH'}
                                       % if HERMES:
                                       %   {'GABAGlx','GSH'}, {'Lac','GSH'}, or {'EtOH','GABA','GSH'}
                                       % if HERCULES:
                                       %   {'GABAGlx','GSH'}
                                       % if phantom data:
                                       %   and MEGA-PRESS: {'GABA'}, {'Glx'}, {'GSH'}, {'Lac'}, or {'EtOH'}
                                       %   and HERMES: {'GABA','GSH'}, {'Glx','GSH'}, {'Lac','GSH'}, or {'EtOH','GABA','GSH'}
    MRS_struct.p.ONOFForder = 'offfirst'; % order of editing pulses; options are 'onfirst' or 'offfirst'
    MRS_struct.p.seqorig = 'JHU'; % origin of Philips patch; options are 'JHU' or 'Philips'
    
% Analysis parameters
    MRS_struct.p.LB = 3; % line-broadening (in Hz)
    MRS_struct.p.water_phase_correction = 1; % 1 = YES; perform eddy current correction on water data
    MRS_struct.p.data_phase_correction = 0; % 1 = YES; perform eddy current correction on metabolite data
    MRS_struct.p.water_removal = 1; % 1 = YES; remove residual water signal using HSVD
    MRS_struct.p.AlignTo = 'RobustSpecReg'; % options are 'RobustSpecReg' (recommended), 'SpecReg', 'SpecRegHERMES',
                                            % 'Cr', 'Cho', 'NAA', 'H2O', 'CrOFF' or 'none'
    MRS_struct.p.Vox = {'vox1','vox2'}; % for naming voxels in PRIAM data, e.g.: {'anterior','posterior'}, {'right','left'}, etc.
    MRS_struct.p.FitResidWater = 0; % 1 = YES, fit the residual water signal in the DIFF spectrum to calculate water suppression factor
    
% Flags
    MRS_struct.p.HERMES   = 0; % 1 = YES, 0 = NO
    MRS_struct.p.HERCULES = 0; % 1 = YES, 0 = NO (if 1, MRS_struct.p.HERMES *must* be set to 1 as well)
    MRS_struct.p.PRIAM    = 0; % 1 = YES, 0 = NO
    MRS_struct.p.phantom  = 0; % 1 = YES (assumes phantom was scanned at room temperature), 0 = NO (for in vivo data)
    MRS_struct.p.mat      = 0; % 1 = YES, save MRS_struct as .mat file
    MRS_struct.p.sdat     = 0; % 1 = YES, save processed difference spectrum as .sdat file (only for Philips SDAT MEGA-PRESS datasets)
    MRS_struct.p.csv      = 0; % 1 = YES, extract useful data from MRS_struct and export to .csv file (applies to GannetFit, GannetSegment and GannetQuantify)
    
end



