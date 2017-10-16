function UTCases = set_up_UT_cases_Perfusion_Phantom;
% run the gt recon

% set_UT_Dir('E')

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {   
    'perfusion\KCL',  'Perfusion_AIF_TwoEchoes_Interleaved_R2_42170_61668703_61668708_324_20170926-154931',      'VD',   'GTPrep_2DT_Perf_AIF_2E_Lin_Mapping_GdPhantom.xml',          'res',          'ref',   'IsmrmrdParameterMap_Siemens_Perfusion.xsl'   ; ...
};

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end
