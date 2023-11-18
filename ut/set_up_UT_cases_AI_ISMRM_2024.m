function UTCases = set_up_UT_cases_AI_ISMRM_2024;
% run the gt recon

% set_UT_Dir('D');

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {

    % --------------------------------------------------------------
    % 
    'neuro', 'meas_MID00094_FID14732_t2_spc_sag_1mm_p2X2', 'VE',   'GTPrep_3DT_Cartesian_AI_FIL.xml',             'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ... 
    'neuro', 'meas_MID00094_FID14732_t2_spc_sag_1mm_p2X2', 'VE',   'Generic_Cartesian_Grappa_SNR.xml',             'res2',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ... 

    '/export/Lab-Kellman/Share/data/LGE/WB/', 'WB_LGE_MOCO_AVE_OnTheFly_41144_01418721_01418731_1929_20230208-164114', 'VE',   'GTPrep_LGE_STCNNT.xml',             'res_ai',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...       
    '/export/Lab-Kellman/Share/data/LGE/WB/', 'WB_LGE_MOCO_AVE_OnTheFly_41144_01418721_01418731_1929_20230208-164114', 'VE',   'GTPrep_LGE.xml',             'res_grappa',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...       


    '/export/Lab-Kellman/Share/data/LGE/WB/', 'WB_LGE_MOCO_AVE_STCNNT_41837_2049151069_2049151078_179_20230904-105942', 'VE',   'GTPrep_LGE_STCNNT.xml',             'res_ai',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...       
    '/export/Lab-Kellman/Share/data/LGE/WB/', 'WB_LGE_MOCO_AVE_STCNNT_41837_2049151069_2049151078_179_20230904-105942', 'VE',   'GTPrep_LGE.xml',             'res_grappa',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...       

    '/export/Lab-Kellman/Share/data/LGE/WB/', 'WB_LGE_MOCO_AVE_STCNNT_41837_2049151069_2049151078_182_20230904-110248', 'VE',   'GTPrep_LGE_STCNNT.xml',             'res_ai',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...       
    '/export/Lab-Kellman/Share/data/LGE/WB/', 'WB_LGE_MOCO_AVE_STCNNT_41837_2049151069_2049151078_182_20230904-110248', 'VE',   'GTPrep_LGE.xml',             'res_grappa',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...       


    'neuro', 'meas_MID00083_FID14721_t1_mprage_1mm_p4_pos50_ACPC_check', 'VE',   'GTPrep_3DT_Cartesian_AI_FIL.xml',             'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ... 

    'LungTSE_rawData', '20181218_124140_meas_MID00542_FID03365_t2_tse_tra_p2_320_trig', 'VE',   'GTPrep_2DT_Cartesian_Lung_TSE.xml',             'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ... 

    'DataForHui/kneeH5/noise', 'noise_20190104_200259_meas_MID00033_FID06439_pd_tse_sag_384', 'VE',   'GTPrep_Measurement_Dependencies_istore.xml',             'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ... 
    
    'DataForHui/spineH5/noise', 'noise_20181207_181828_meas_MID00150_FID01036_t2_tse_sag_p2', 'VE',   'GTPrep_Measurement_Dependencies_istore.xml',             'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ... 
    'DataForHui/spineH5/imagedata', '20181207_181828_meas_MID00150_FID01036_t2_tse_sag_p2', 'VE',   'Generic_Cartesian_Grappa_SNR_STCNNT.xml',             'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ... 

    % ------------------------------------------------------------
};

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end
