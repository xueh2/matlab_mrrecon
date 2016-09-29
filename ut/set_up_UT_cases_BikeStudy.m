function UTCases = set_up_UT_cases_BikeStudy;
% run the gt recon

set_UT_Dir('D')

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = { 
            'cine\BikeStudy\CNMC2359-VO813',                 '20150826_16h11m36s_59731',          'VD',   'Generic_Cartesian_NonLinear_Spirit_RealTimeCine_Cloud.xml',  'nl_spirit_cloud_res_db1_noise_floor_0p0015', 'nl_spirit_cloud_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...                     
            'cine\BikeStudy\CNMC2359-VO813',                 '20150826_16h36m28s_59749',          'VD',   'GTPrep_2DT_RealTimeCine_PF_Handling_PhysioInterp.xml',  'grappa', 'grappa_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...                                 
            
            'cine\BikeStudy\CNMC2359-VO817',                 '20150828_14h19m02s_60279',          'VD',   'Generic_Cartesian_NonLinear_Spirit_RealTimeCine_Cloud.xml',  'nl_spirit_cloud_res_db1_noise_floor_0p0015', 'nl_spirit_cloud_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...                                 
            'cine\BikeStudy\CNMC2359-VO817',                 '20150828_14h32m35s_60288',          'VD',   'CMR_2DT_RTCine_KspaceBinning_Cloud.xml',  'binning_cloud_res', 'binning_cloud_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...                                 
            
            'cine\BikeStudy\CNMC2359-VO933',                 '20160218_17h49m39s_2045',          'VD',   'GTPrep_2DT_RealTimeCine_PF_Handling_PhysioInterp.xml',  'grappa', 'grappa_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...                                 
            
            'cine\BikeStudy\VO583',                 '20151210_15h16m50s_72947',          'VD',   'Generic_Cartesian_NonLinear_Spirit_RealTimeCine_Cloud.xml',  'nl_spirit_cloud_res_db1_noise_floor_0p0015', 'nl_spirit_cloud_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...                     
            'cine\BikeStudy\VO583',                 '20151210_15h16m00s_72945',          'VD',   'GTPrep_2DT_RealTimeCine_PF_Handling_PhysioInterp.xml',  'grappa', 'grappa_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...                                 
            'cine\BikeStudy\VO583',                 '20151210_15h13m54s_72944',          'VD',   'CMR_2DT_RTCine_KspaceBinning_Cloud.xml',  'grappa', 'grappa_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...                                 
            
    };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end
