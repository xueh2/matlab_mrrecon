function UTCases = set_up_UT_cases_Cloud;
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = { % test
            'binning',                  'meas_MID836_PK_rt_test_1slice_FID22517',              'VD',   'CMR_2DT_RTCine_KspaceBinning_Cloud.xml',         'binning_cloud_res',              'binning_cloud_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...    
            'binning',                  'meas_MID838_PK_rt_test_2slice_FID22519',              'VD',   'CMR_2DT_RTCine_KspaceBinning_Cloud.xml',         'binning_cloud_res',              'binning_cloud_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...    
            
            'binning',                  'Cmr_Binning_Cloud_42110_10890773_10890782_2694_20201009-171402',              'VD',   'CMR_2DT_RTCine_KspaceBinning_Cloud.xml',         'binning_cloud_res',              'binning_cloud_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...    
            
            % debug, 20160622
            'binning',                  'meas_MID00060_FID05129_SAX_RT_Binning_8sec_Lighthouse_(Baseline)', 'VD',   'CMR_2DT_RTCine_KspaceBinning_Cloud.xml',         'binning_cloud_res',              'binning_cloud_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...    
            
            % debug, 20160715
            'binning\CNMC_20160715',   'meas_MID00062_FID23728_SAX_R3_RT_Binning_256_CLOUD_9sl',              'VD',   'CMR_2DT_RTCine_KspaceBinning_Cloud.xml',                                  'binning_cloud_res',            'binning_cloud_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...   
            
            'binning',                  'Cmr_Binning_Cloud_41548_220125221_220125228_42_20210319-085028', 'VD',   'CMR_2DT_RTCine_KspaceBinning_Cloud.xml',         'binning_cloud_res',              'binning_cloud_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...    
    };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end
