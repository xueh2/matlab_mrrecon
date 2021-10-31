function UTCases = set_up_UT_cases_PreScan;
% UTCases = set_up_UT_cases_PreScan;

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {
    'PreScan',               'ISMRMRD_Noise_dependency_41837_4052642_4052651_974_20210305-121822',            'VD',   'default_measurement_dependencies_Noise_CoilSen_SCC.xml',         'res',              'ref',   'IsmrmrdParameterMap_Siemens_SashaT1T2_HC.xsl'   ; ...    
    'PreScan',               'RT_Cine_LIN_41837_4052642_4052651_990_20210305-122857',                         'VD',   'GTPrep_2DT_Cartesian_SCC.xml',         'res',              'ref',   'IsmrmrdParameterMap_Siemens_SashaT1T2_HC.xsl'   ; ...    
    }

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end
