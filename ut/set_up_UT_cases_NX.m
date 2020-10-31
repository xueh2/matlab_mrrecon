function UTCases = set_up_UT_cases_NX;
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {

    'VIDA_TEST_DATA',         'meas_MID00174_FID04926_de_tfl_high_res_psir_LG',         'VD',   'GTPrep_2DT_LGE_MOCO_AVE_Denoise_OnTheFly.xml',                          'moco_ave_res_fil',            'moco_ave_ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
    
    
    'NXVA11A',         'meas_MID00021_FID29399_BEAT_epi_DS_AIF',         'NX',   'GTPrep_2DT_Perf_AIFR3_2E_Lin_Mapping_MBF_MBV_Mask_CMR_View_OFFLINE.xml',                          'moco_ave_res_fil',            'moco_ave_ref',   'IsmrmrdParameterMap_Siemens_Perfusion_NX.xsl'   ; ...
    'NXVA11A',         'meas_MID00025_FID29299_BEAT_epi_DS_AIF',         'NX',   'GTPrep_2DT_Perf_AIFR3_2E_Lin_Mapping_MBF_MBV_Mask_CMR_View_OFFLINE.xml',                          'moco_ave_res_fil',            'moco_ave_ref',   'IsmrmrdParameterMap_Siemens_Perfusion_NX.xsl'   ; ...
    
    'NXVA11A',         'meas_MID00031_FID30684_PSIR_MOCO_8_avg_9_slice_GT',         'NX',   'GTPrep_2DT_LGE_MOCO_AVE_OnTheFly.xml',                          'moco_ave_res_fil',            'moco_ave_ref',   'IsmrmrdParameterMap_Siemens_NX.xsl'   ; ...

    'NormalOrien/BARTS_Perfusion_20191011',         'meas_MID00314_FID259832_SSFP_Test_Perf_MBF_AI',         'NX',   'GTPrep_2DT_Perf_AIFR3_2E_Lin_Mapping_MBF_MBV_Mask_CMR_View_OFFLINE.xml',                          'moco_ave_res_fil',            'moco_ave_ref',   'IsmrmrdParameterMap_Siemens_Perfusion_NX.xsl'   ; ...
    
    'NXVA11A',         'Perfusion_AIFR3_2E_Interleaved_175749_00002607024.2607061_00002607024.2607061_322_20191022-125558',         'NX',   'GTPrep_2DT_Perf_AIFR3_2E_Lin_Mapping_MBF_MBV_Mask_CMR_View_OFFLINE.xml',                          'res',            'ref',   'IsmrmrdParameterMap_Siemens_Perfusion_NX.xsl'   ; ...
    
    'D:\data\test\XA20',         'meas_MID00032_FID00510_SSFP_BEAT_epi_DS_AIF',         'NX',   'GTPrep_2DT_Perf_AIF_2E_Lin_Mapping_MBF_MBV_Mask_AI_CMR_View_OFFLINE.xml',                          'res',            'ref',   'IsmrmrdParameterMap_Siemens_Perfusion_NX20A.xsl'   ; ...
    
    };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end
