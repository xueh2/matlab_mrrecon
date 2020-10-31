function UTCases = set_up_UT_cases_BinningCine;
% run the gt recon

% set_UT_Dir('D')

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {  
    
    'binning',                  'meas_MID836_PK_rt_test_1slice_FID22517',              'VD',   'GTPrep_2DT_RTCine_KspaceBinning.xml',      'grappa_res',               'grappa_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...    
    'binning',                  'meas_MID836_PK_rt_test_1slice_FID22517',              'VD',   'CMR_2DT_RTCine_KspaceBinning.xml',         'generic_res',              'generic_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...    
    
    'binning',                  'meas_MID838_PK_rt_test_2slice_FID22519',              'VD',   'CMR_2DT_RTCine_KspaceBinning.xml',         'generic_res',              'generic_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...    
    
    % 2CH
    'binning',                  'Cmr_Binning_Cloud_42537_96833260_96833270_463_20200910-141550',     'VD',   'CMR_2DT_Binning_Fil_LAX_GLS_AI.xml',         'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...        
    
    'binning',                  'Cmr_Binning_Cloud_42537_96833260_96833270_463_20200910-141550',     'VD',   'CMR_2DT_Binning_Fil_LAX_GLS_AI.xml',         'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...        
    'binning',                  'Cmr_Binning_Cloud_41672_33673004_33673012_845_20201009-141458',     'VD',   'CMR_2DT_Binning_Fil_LAX_GLS_AI.xml',         'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...        
    'binning',                  'Cmr_Binning_Cloud_42110_10890773_10890782_2694_20201009-171402',     'VD',   'CMR_2DT_Binning_Fil_LAX_GLS_AI.xml',         'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...        
    'binning',                  'Cmr_Binning_Cloud_42110_10890773_10890782_2696_20201009-171435',     'VD',   'CMR_2DT_Binning_Fil_LAX_GLS_AI.xml',         'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...        
    'binning',                  'Cmr_Binning_Cloud_42110_10890773_10890782_2698_20201009-171501',     'VD',   'CMR_2DT_Binning_Fil_LAX_GLS_AI.xml',         'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...        
   };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end
