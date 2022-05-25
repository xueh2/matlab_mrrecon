function UTCases = set_up_UT_cases_PreScan;
% UTCases = set_up_UT_cases_PreScan;

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {
    'PreScan',               'ISMRMRD_Noise_dependency_41837_4052642_4052651_974_20210305-121822',            'VD',   'default_measurement_dependencies_Noise_CoilSen_SCC.xml',         'res',              'ref',   'IsmrmrdParameterMap_Siemens_SashaT1T2_HC.xsl'   ; ...    
    'PreScan',               'RT_Cine_LIN_41837_4052642_4052651_990_20210305-122857',                         'VD',   'Generic_RTCine_SCC.xml',         'res',              'ref',   'IsmrmrdParameterMap_Siemens_SashaT1T2_HC.xsl'   ; ...    

    'PreScan/meas',    '20220225_113745_meas_MID00026_FID38202_Ch2_Localizer_simple_PAT2_SCC',            'VD',   'default_measurement_dependencies_Noise_CoilSen_SCC.xml',         'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...    
    'PreScan/meas',    '20220225_113745_meas_MID00026_FID38202_Ch2_Localizer_simple_PAT2',            'VD',   'Generic_RTCine_SCC.xml',         'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...    

    'PreScan/meas',    '20220225_113854_meas_MID00035_FID38211_Ch4_Localizer_simple_PAT2_SCC',            'VD',   'default_measurement_dependencies_Noise_CoilSen_SCC.xml',         'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...    
    'PreScan/meas',    '20220225_113854_meas_MID00035_FID38211_Ch4_Localizer_simple_PAT2',            'VD',   'Generic_RTCine_SCC.xml',         'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...    

    %'PreScan/meas',    '20220225_113854_meas_MID00035_FID38211_Ch4_Localizer_simple_PAT2',            'VD',   'default_measurement_dependencies_Noise_CoilSen_SCC.xml',         'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...    
    }

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end

command = 'gadgetron_ismrmrd_client -f T:\Lab-Kellman\Share\data\PreScan\meas\20220225_113745_meas_MID00026_FID38202_Ch2_Localizer_simple_PAT2_SCC\20220225_113745_meas_MID00026_FID38202_Ch2_Localizer_simple_PAT2_SCC.h5 -C D:\gtuser\mrprogs\install_master_debug\share\gadgetron\config\default_measurement_dependencies_Noise_CoilSen_SCC.xml -a localhost -p 9002 -F hdr -G default_measurement_dependencies_Noise_CoilSen_SCC.xml -o ref_20220408_135432.h5'

command_im = 'gadgetron_ismrmrd_client -f T:\Lab-Kellman\Share\data\PreScan\meas\20220225_113745_meas_MID00026_FID38202_Ch2_Localizer_simple_PAT2\20220225_113745_meas_MID00026_FID38202_Ch2_Localizer_simple_PAT2.h5 -C D:\gtuser\mrprogs\install_master_debug\share\gadgetron\config\Generic_RTCine_SCC.xml -a localhost -p 9002 -F hdr -G Generic_RTCine_SCC.xml -o ref_20220408_135432.h5'

