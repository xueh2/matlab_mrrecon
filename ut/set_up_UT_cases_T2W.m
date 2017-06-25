function UTCases = set_up_UT_cases_T2W;
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {
    'T2W',         'meas_MID00057_FID25460_QQQ_4ch_CV_T2w_567C_MOCO',         'VD',   'GTPrep_2DT_T2W_MOCO_AVE.xml',                                         'moco_ave_res',            'moco_ave_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
    'T2W',         'meas_MID00123_FID43785_gt_CV_T2w_segmented',         'VD',   'GTPrep_2DT_T2W_MOCO_AVE_OnTheFly.xml',                                         'moco_ave_res',            'moco_ave_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
    'T2W',         'meas_MID00124_FID43786_gt_T2w_1sl_8ave_MOCO',         'VD',   'GTPrep_2DT_T2W_MOCO_AVE_OnTheFly.xml',                                         'moco_ave_res',            'moco_ave_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
    'T2W',         'meas_MID00371_FID108006_SAX_gt_T2w_9sl_8ave_MOCO',         'VD',   'GTPrep_2DT_T2W_MOCO_AVE_OnTheFly.xml',                                         'moco_ave_res',            'moco_ave_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
    
    't2w\T2pSSFP_Artifact',         'meas_MID00309_FID64148_FB_MOCO_T2SSFP_1sl_8av',         'VD',   'GTPrep_2DT_T2W_MOCO_AVE_OnTheFly.xml',                                         'moco_ave_res',            'moco_ave_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...

    };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];

end
