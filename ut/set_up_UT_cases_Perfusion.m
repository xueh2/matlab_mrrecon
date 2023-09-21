function UTCases = set_up_UT_cases_Perfusion;
% run the gt recon

% set_UT_Dir('E')

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {   

    '/export/Lab-Kellman/Share/data/perfusion/UVA/',  'meas_MID00223_FID196536_GRE_qPerf_MBF_AI_STRESS',  'NX20',   'GTPrep_2DT_Perf_AIF_2E_AI_NR_STCNNT_OFFLINE.xml',          'res',          'grappa',   'wip_071_qPerf_IsmrmrdParameterMap_Siemens_Perfusion_NX20.xsl'   ; ...                            
    '/export/Lab-Kellman/Share/data/perfusion/UVA/',  'meas_MID00225_FID196538_GRE_qPerf_MBF_AI_REST',  'NX20',   'GTPrep_2DT_Perf_AIF_2E_AI_NR_STCNNT_OFFLINE.xml',          'res',          'grappa',   'wip_071_qPerf_IsmrmrdParameterMap_Siemens_Perfusion_NX20.xsl'   ; ...                            

    
    };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end
