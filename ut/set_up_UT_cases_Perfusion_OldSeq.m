function UTCases = set_up_UT_cases_Perfusion_OldSeq;
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {
    'FLASH_Perfusion_HighRes',  '20120809_11h55m25s_8477',      'VD',   'GTPrep_2DT_Perfusion_AIF.xml',          'grappa_res',          'grappa_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
    };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end
