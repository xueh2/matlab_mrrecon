function UTCases = set_up_UT_cases_DicomStorage;
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {
    'LGE\SCC_Test', '20150302_10h50m25s_66537_h5', 'VD',   'GTPrep_2DT_LGE_MOCO_AVE_OnTheFly_ismrmrd_storage.xml',   'dicom_res', 'dicom_ref', 'IsmrmrdParameterMap_Siemens.xsl'   ; ...       
    'GT_UT\CoilQA', 'meas_MID00202_FID109023_Coil_QA', 'VD',   'GTPrep_2DT_Cartesian_Dicom.xml',   'dicom_res', 'dicom_ref', 'IsmrmrdParameterMap_Siemens.xsl'   ; ...       
   };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end
