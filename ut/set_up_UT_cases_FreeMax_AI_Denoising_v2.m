function UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;
% run the gt recon

% set_UT_Dir('D');

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {

    '/export/Lab-Xue/data/mri_raw_data/freemax/20230630_NV_AI',         'meas_MID00155_FID07561_G25_4CH_CINE_256_R2',              'NX',   'GTPrep_2DT_RetroCine_STCNNT_offline.xml',             'res',              'ref',   'wip_070_fire_IsmrmrdParameterMap_Siemens.xsl'   ; ...    
    '/export/Lab-Xue/data/mri_raw_data/freemax/20230630_NV_AI',         'meas_MID00156_FID07562_G25_2CH_CINE_256_R2',              'NX',   'GTPrep_2DT_RetroCine_STCNNT_offline.xml',             'res',              'ref',   'wip_070_fire_IsmrmrdParameterMap_Siemens.xsl'   ; ...    

    '/export/Lab-Xue/data/mri_raw_data/freemax/20230630_NV_AI',         'meas_MID00159_FID07565_G25_2CH_CINE_256_R3',              'NX',   'GTPrep_2DT_RetroCine_STCNNT_offline.xml',             'res',              'ref',   'wip_070_fire_IsmrmrdParameterMap_Siemens.xsl'   ; ...    

    '/export/Lab-Xue/data/mri_raw_data/freemax/20230630_NV_AI',         'meas_MID00163_FID07569_G25_4CH_CINE_256_R4',              'NX',   'GTPrep_2DT_RetroCine_STCNNT_offline.xml',             'res',              'ref',   'wip_070_fire_IsmrmrdParameterMap_Siemens.xsl'   ; ...       
    '/export/Lab-Xue/data/mri_raw_data/freemax/20230630_NV_AI',         'meas_MID00164_FID07570_G25_2CH_CINE_256_R4',              'NX',   'GTPrep_2DT_RetroCine_STCNNT_offline.xml',             'res',              'ref',   'wip_070_fire_IsmrmrdParameterMap_Siemens.xsl'   ; ...    
    
    '/export/Lab-Xue/data/mri_raw_data/freemax/20230630_NV_AI',         'meas_MID00176_FID07582_G25_Perfusion_trufi_sr_tpat_4_256res',              'NX',   'GTPrep_2DT_Perf_AIF_2E_AI_NR_STCNNT.xml',             'res',              'ref',   'wip_070_fire_IsmrmrdParameterMap_Siemens.xsl'   ; ...    
    '/export/Lab-Xue/data/mri_raw_data/freemax/20230630_NV_AI',         'meas_MID00175_FID07581_G25_Perfusion_trufi_sr_tpat_3_192res',              'NX',   'GTPrep_2DT_Perf_AIF_2E_AI_NR_STCNNT.xml',             'res',              'ref',   'wip_070_fire_IsmrmrdParameterMap_Siemens.xsl'   ; ...    
    };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end
