function UTCases = set_up_UT_cases_PerfusionCloud;
% run the gt recon

% set_UT_Dir('E')

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {   
    'perfusion/cloud',  'Perfusion_AIF_2E_NL_Cloud_42170_072714971_072714980_550_20180718-175707',      'VD',   'GTPrep_2DT_Perf_AIF_2E_Lin_Mapping_MBF_MBV_Mask_Gateway_istore.xml',          'cloud_grappa_flow_res',          'cloud_flow_ref',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ; ...
    'perfusion/cloud',  'Perfusion_AIF_2E_NL_Cloud_42170_072714971_072714980_550_20180718-175707',      'VD',   'GTPrep_2DT_Perf_AIF_2E_NL_Mapping_MBF_MBV_Mask_Gateway.xml',          'cloud_flow_res',          'cloud_flow_ref',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ; ...
    
};

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end
