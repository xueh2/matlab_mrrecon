function UTCases = set_up_UT_cases_Perfusion;
% run the gt recon

% set_UT_Dir('E')

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {   

    'perfusion/Barts',  'Perfusion_AIF_TwoEchoes_Interleaved_R2_42110_179930079_179930088_2790_20220202-172819',  'VE',   'GTPrep_2DT_Perf_AIF_2E_AI_NR_STCNNT_OFFLINE.xml',          'res',          'grappa',   'wip_071_qPerf_IsmrmrdParameterMap_Siemens_Perfusion_NX20.xsl'   ; ...                            
    'perfusion/Barts',  'Perfusion_AIF_TwoEchoes_Interleaved_R2_42110_179930079_179930088_2793_20220202-173329',  'VE',   'GTPrep_2DT_Perf_AIF_2E_AI_NR_STCNNT_OFFLINE.xml',          'res',          'grappa',   'wip_071_qPerf_IsmrmrdParameterMap_Siemens_Perfusion_NX20.xsl'   ; ...                            
    'perfusion/Barts',  'Perfusion_AIF_TwoEchoes_Interleaved_R2_66016_1340151873_1340151882_95_20201001-102249',  'VE',   'GTPrep_2DT_Perf_AIF_2E_AI_NR_STCNNT_OFFLINE.xml',          'res',          'grappa',   'wip_071_qPerf_IsmrmrdParameterMap_Siemens_Perfusion_NX20.xsl'   ; ...                            
    'perfusion/Barts',  'Perfusion_AIF_TwoEchoes_Interleaved_R2_66016_1340151873_1340151882_98_20201001-102934',  'VE',   'GTPrep_2DT_Perf_AIF_2E_AI_NR_STCNNT_OFFLINE.xml',          'res',          'grappa',   'wip_071_qPerf_IsmrmrdParameterMap_Siemens_Perfusion_NX20.xsl'   ; ...                            

    };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end
