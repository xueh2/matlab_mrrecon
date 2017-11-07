function UTCases = set_up_UT_cases_BART;
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {
    'Barts', 'meas_MID00168_FID39823_R3_RT_Cine_res160_46ms_NL_BART_Proto1',  'VD',   'BART_Recon.xml',  'grappa_res',  'grappa_ref',   'IsmrmrdParameterMap_Siemens.xsl'  ; ...   
    'Barts', 'meas_MID00168_FID39823_R3_RT_Cine_res160_46ms_NL_BART_Proto1',  'VD',   'BART_Recon_Cloud.xml',  'cloud_res',  'cloud_ref',   'IsmrmrdParameterMap_Siemens.xsl'  ; ...   
    'cine', 'meas_MID00165_FID39820_R3_RT_Cine_Linear_Res160_35ms',  'VD',   'BART_Recon_Cloud.xml',  'cloud_res',  'cloud_ref',   'IsmrmrdParameterMap_Siemens.xsl'  ; ...   
    'cine', 'meas_MID00042_FID96308_cine_SAX_realtime_gt_TPAT4',  'VD',   'BART_Recon_Cloud.xml',  'cloud_res',  'cloud_ref',   'IsmrmrdParameterMap_Siemens.xsl'  ; ...   
   };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end
