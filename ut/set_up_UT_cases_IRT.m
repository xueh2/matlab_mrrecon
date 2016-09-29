function UTCases = set_up_UT_cases_IRT
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {
'GT_UT/IRT',  'meas_MID00137_FID24285_2015_12_7_BEAT_IRTTT_1030_GT', 'VD',   'grappa_device.xml',      'res',           'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...   
'GT_UT/IRT',  'meas_MID00137_FID24285_2015_12_7_BEAT_IRTTT_1030_GT', 'VD',   'grappa_device_cpu.xml',  'res_cpu',       'ref_cpu',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...   
};

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end
