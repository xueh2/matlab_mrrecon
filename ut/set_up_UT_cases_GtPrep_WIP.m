function UTCases = set_up_UT_cases_GtPrep_WIP;
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {
    

    'gtprep_data/gtprep/Perfusion/wips',       'Perfusion_AIF_TwoEchoes_Interleaved_R2_42110_3662928_3662937_656_20210127-120816',  'VD',   'Gadgetron_QPerf_WIP_AI.xml',           'res',           'ref',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ; ...
    'gtprep_data/gtprep/Perfusion/wips',       'Perfusion_AIF_TwoEchoes_Interleaved_R2_42110_3662928_3662937_659_20210127-121722',  'VD',   'Gadgetron_QPerf_WIP_AI.xml',           'res',           'ref',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ; ...
    };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

disp('set_up_UT_cases_GtPrep_WIP')

if(nargout<=0)
    UTCases = [];
end
