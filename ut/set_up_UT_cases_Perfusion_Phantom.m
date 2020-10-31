function UTCases = set_up_UT_cases_Perfusion_Phantom;
% run the gt recon

% set_UT_Dir('E')

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {   
    'perfusion/KCL',  'Perfusion_AIF_TwoEchoes_Interleaved_R2_42170_61668703_61668708_324_20170926-154931',      'VD',   'GTPrep_2DT_Perf_AIF_2E_Lin_Mapping_GdPhantom.xml',          'res',          'ref',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ; ...
    
    'perfusion/Perf_no_AIF',  'meas_MID01771_FID134128_SSFP_Perf_no_map',      'VD',   'GTPrep_2DT_Perf_NoAIF_Mapping.xml',          'res',          'ref',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ; ...
    'perfusion/Perf_no_AIF',  'meas_MID01772_FID134129_SSFP_Perf_no_map',      'VD',   'GTPrep_2DT_Perf_NoAIF_Mapping.xml',          'res',          'ref',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ; ...
    
    'perfusion',  'meas_MID00293_FID42018_SSFP_PERFUSION_GADGETRON_AIF2eTPAT_1slc_96res',      'VD',   'GTPrep_2DT_Perf_AIFR3_2E_Lin_Mapping_MBF_MBV_Mask_AI_CMR_View_OFFLINE.xml',          'res',          'ref',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ; ...
};

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end
