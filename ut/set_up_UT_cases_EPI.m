function UTCases = set_up_UT_cases_EPI;
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {
    'epi', 'meas_MID00087_FID10658_amri_ep2d_bold_Dante_ish_72x54x50_3000msTR_28msTE',  'VB',   'epi.xml',  'grappa_res',  'grappa_ref',   'IsmrmrdParameterMap_Siemens_EPI.xsl'  ; ...   
    'epi', 'meas_MID00087_FID10658_amri_ep2d_bold_Dante_ish_72x54x50_3000msTR_28msTE',  'VB',   'Generic_Cartesian_Grappa_EPI.xml',  'grappa_res',  'grappa_ref',   'IsmrmrdParameterMap_Siemens_EPI.xsl'  ; ...   
    'epi', 'meas_MID00063_FID10883_amri_ep2d_bold_Dante_ish_96x72x50R2_3000msTR_24msTE',  'VD',   'Generic_Cartesian_Grappa_EPI_AVE.xml',  'grappa_res',  'grappa_ref',   'IsmrmrdParameterMap_Siemens_EPI.xsl'  ; ...   
    'epi', 'meas_MID00065_FID10885_amri_ep2d_bold_Dante_ish_96x72x50R2_3000msTR_24msTE_R2',  'VD',   'Generic_Cartesian_Grappa_EPI_AVE.xml',  'grappa_res',  'grappa_ref',   'IsmrmrdParameterMap_Siemens_EPI_FLASHREF.xsl'  ; ...   
    'epi', 'meas_MID00275_FID11082_amri_ep2d_bold_Dante_ish_96x72x50_Gadgetron_test_R1_WRONG_TR',  'VD',   'Generic_Cartesian_Grappa_EPI_AVE.xml',  'grappa_res',  'grappa_ref',   'IsmrmrdParameterMap_Siemens_EPI_FLASHREF.xsl'  ; ...   
    'epi', 'meas_MID01349_FID12150_amri_ep2d_bold_96x72x5_R2_16avg_gadgetron',  'VD',   'Generic_Cartesian_Grappa_EPI_AVE.xml',  'grappa_res',  'grappa_ref',   'IsmrmrdParameterMap_Siemens_EPI_FLASHREF.xsl'  ; ...   
   };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end
