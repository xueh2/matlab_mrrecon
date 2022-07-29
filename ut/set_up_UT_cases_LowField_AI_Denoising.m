function UTCases = set_up_UT_cases_LowField_AI_Denoising;
% run the gt recon

% set_UT_Dir('D');

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {
    
    'FreeMax/20220718_NV_cardiac',         'CINE_PAT2_1slice_RR_206026_30000022071809030959400000008_30000022071809030959400000008_3000018_20220729-132033',              'NX',   'GTPrep_2DT_RetroCine_CNNT_SCC_GLS_Seg_AI_istore.xml',             'res',              'ref',   'IsmrmrdParameterMap_Siemens_Perfusion_NX.xsl'   ; ...

    'FreeMax/20220718_NV_cardiac',         'meas_MID00037_FID05411_CINE_PAT3_1slice',              'NX',   'Generic_RTCine_CNNT.xml',             'res',              'ref',   'IsmrmrdParameterMap_Siemens_Perfusion_NX.xsl'   ; ...
   
    'LowField/CNNT/binning',               '20200129_074228_meas_MID00100_FID60958_4CH_BinningCine_R3_192res_cloud',              'VD',   'Generic_RTCine_CNNT.xml',             'res',              'ref',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
    'LowField/CNNT/binning',               '20200129_074303_meas_MID00101_FID60959_3CH_BinningCine_R3_192res_cloud',              'VD',   'Generic_RTCine_CNNT.xml',             'res',              'ref_fil',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...    
    'LowField/CNNT/binning',               '20200129_074344_meas_MID00102_FID60960_2CH_BinningCine_R3_192res_cloud',              'VD',   'Generic_RTCine_CNNT.xml',             'res',              'ref_fil',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...

    'LowField/CNNT/binning',               '20190701_102821_meas_MID00037_FID24106_4Ch_BinningCine_R3_256res_BW751_35s',              'VD',   'Generic_RTCine_CNNT.xml',             'res',              'ref_fil',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
    'LowField/CNNT/binning',               '20190701_102923_meas_MID00038_FID24107_2Ch_BinningCine_R3_256res_BW751_35s',              'VD',   'Generic_RTCine_CNNT.xml',             'res',              'ref_fil',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
    'LowField/CNNT/binning',               '20190701_102237_meas_MID00036_FID24105_SAX_BinningCine_R3_256res_BW751_35s',              'VD',   'Generic_RTCine_CNNT.xml',             'res',              'ref_fil',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...

    'LowField/CNNT/neuro',               '20190806_190618_meas_MID00129_FID32022_t1_se_tra_manual_denoising',                   'VD',   'GTPrep_2DT_Cartesian_GFactor_Offline.xml',             'res',              'ref_fil',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
    'LowField/CNNT/neuro',               '20190809_190237_meas_MID00043_FID32676_Post_t1_se_tra_manual_denoising',              'VD',   'GTPrep_2DT_Cartesian_GFactor_Offline.xml',             'res',              'ref_fil',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...

    'LowField/CNNT/spine',               '20200305_212932_meas_MID00031_FID67000_t2_tse_sag_c-spine',              'VD',   'GTPrep_2DT_AVE_CNNT.xml',             'res',              'ref_fil',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
    'LowField/CNNT/spine',               '20200305_213517_meas_MID00037_FID67006_t2_tse_sag_t-spine',              'VD',   'GTPrep_2DT_AVE_CNNT.xml',             'res',              'ref_fil',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
    'LowField/CNNT/spine',               '20200305_214034_meas_MID00040_FID67009_t2_tse_sag_l-spine',              'VD',   'GTPrep_2DT_AVE_CNNT.xml',             'res',              'ref_fil',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
    'LowField/CNNT/spine',               '20190225_213935_meas_MID00043_FID00868_t2_tse_sag_c-spine',              'VD',   'GTPrep_2DT_AVE_CNNT.xml',             'res',              'ref_fil',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...   
    'LowField/CNNT/spine',               '20181207_213119_meas_MID00183_FID01069_t2_tse_sag_c-spine',              'VD',   'GTPrep_2DT_AVE_CNNT.xml',             'res',              'ref_fil',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...
%     
%       
%     20181207_213657_meas_MID00186_FID01072_t2_tse_sag_t-spine  
%     20181207_220923_meas_MID00202_FID01088_t2_tse_sag_l-spine  
%     20190225_213935_meas_MID00043_FID00868_t2_tse_sag_c-spine  
%     20190314_170101_meas_MID00041_FID04242_t2_tse_sag_l-spine  
%     20190611_172307_meas_MID00041_FID19145_t2_tse_sag_c-spine  
%     20190611_173611_meas_MID00064_FID19168_t2_tse_sag_t-spine  
%     20190611_174100_meas_MID00070_FID19174_t2_tse_sag_l-spine	
%     20190627_163244_meas_MID00056_FID23471_t2_tse_sag_t-spine  
%     20190627_165614_meas_MID00080_FID23495_t2_tse_sag_l-spine  
%     20190719_221338_meas_MID00104_FID29096_t2_tse_sag_c-spine  
%     20200121_212439_meas_MID00034_FID59793_t2_tse_sag_c-spine  
%     20200121_213058_meas_MID00046_FID59805_t2_tse_sag_t-spine  
%     20200121_213620_meas_MID00052_FID59811_t2_tse_sag_l-spine  
%     20200211_141456_meas_MID00022_FID62939_t2_tse_sag_c-spine  
%     20200305_212932_meas_MID00031_FID67000_t2_tse_sag_c-spine	
%     20200305_213517_meas_MID00037_FID67006_t2_tse_sag_t-spine  
%     20200305_214034_meas_MID00040_FID67009_t2_tse_sag_l-spine
    
    };

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end
