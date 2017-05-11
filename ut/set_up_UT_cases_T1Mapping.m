function UTCases = set_up_UT_cases_T1Mapping;
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {
    'sasha\SASHA_PHANTOM', 'meas_MID00613_FID40909_pre_SASHA_10pt_TS600_256', 'VD', 'CMR_2DT_T1Mapping_SASHA.xml', 'grappa_res', 'grappa_ref', 'IsmrmrdParameterMap_Siemens_T1Mapping_SASHA.xsl'   ; ...   
    'sasha\SASHA_PHANTOM', 'meas_MID00614_FID40910_pre_SASHA_10pt_variable_256', 'VD', 'CMR_2DT_T1Mapping_SASHA.xml', 'grappa_res', 'grappa_ref', 'IsmrmrdParameterMap_Siemens_T1Mapping_SASHA.xsl'   ; ...   
    
    'sasha', 'meas_MID01026_FID67745_pre_SASHA_2_1_0_8_1041C_V30_mrtime', 'VD', 'CMR_2DT_T1Mapping_SASHA.xml', 'grappa_res', 'grappa_ref', 'IsmrmrdParameterMap_Siemens_T1Mapping_SASHA.xsl'   ; ...   
    };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];

end
