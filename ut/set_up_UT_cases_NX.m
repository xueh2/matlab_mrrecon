function UTCases = set_up_UT_cases_NX;
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {

    'VIDA_TEST_DATA',         'meas_MID00174_FID04926_de_tfl_high_res_psir_LG',         'VD',   'GTPrep_2DT_LGE_MOCO_AVE_Denoise_OnTheFly.xml',                          'moco_ave_res_fil',            'moco_ave_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
    };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end
