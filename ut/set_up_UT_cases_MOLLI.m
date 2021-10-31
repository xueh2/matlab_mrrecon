function UTCases = set_up_UT_cases_MOLLI;
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {
     
    
    'molli',                          '20100330_10h33m11s_5562',                                  'VD', 'CMR_2DT_T1Mapping_MOLLI.xml', 'res',       'ref',   'IsmrmrdParameterMap_Siemens_T1Mapping_MOLLI.xsl'   ; ...    
    'molli',                          '20100330_10h33m11s_5562',                                  'VD', 'GTPrep_2DT_MOLLI_Offline.xml', 'res',       'ref',   'IsmrmrdParameterMap_Siemens_T1Mapping_MOLLI.xsl'   ; ...    
    'molli',                          'meas_MID00763_FID130687_pre_MOLLI_5s3s3s_256',             'VD', 'MOLLI_T1_Moco.xml', 'res',       'ref',   'IsmrmrdParameterMap_Siemens_T1Mapping_MOLLI.xsl'   ; ...    
    
    'molli',                          'post_MOLLI_4s(1s)3s(1s)2s_256_42363_245912650_245912655_58_20210223-095038',             'VD', 'MOLLI_T1_Moco_istore.xml', 'res',       'ref',   'IsmrmrdParameterMap_Siemens_T1Mapping_MOLLI.xsl'   ; ...    
   };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end


