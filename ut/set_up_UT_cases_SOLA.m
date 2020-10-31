function UTCases = set_up_UT_cases_SOLA;
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {    
    '20200313_LUND_SOLA_XA11B_SNR_TESTDATA/meas_MID00185_FID02098_localizer_128_noacc',           'meas_MID00185_FID02098_localizer_128_noacc',  'NX',   'Generic_Cartesian_Grappa_SNR.xml', 'res',       'ref',   'IsmrmrdParameterMap_Siemens_NX.xsl'   ; ...        
    
    '20200313_LUND_SOLA_XA11B_SNR_TESTDATA/meas_MID00184_FID02097_localizer_128_PAT_4_GREseparate',           'meas_MID00184_FID02097_localizer_128_PAT_4_GREseparate',  'NX',   'Generic_Cartesian_Grappa_SNR.xml', 'res',       'ref',   'IsmrmrdParameterMap_Siemens_NX.xsl'   ; ...            
    
    'SOLA\LUND\20200312',           'meas_MID00250_FID01559_stressBEAT_epi_DS_AIF',  'NX',   'GTPrep_2DT_Perf_AIFR3_2E_Lin_Mapping_MBF_MBV_Mask_AI_CMR_View_OFFLINE.xml', 'res',       'ref',   'IsmrmrdParameterMap_Siemens_Perfusion_NX.xsl'   ; ...            
        
    'SOLA/20201001_Marecav_pat1',           'meas_MID00571_FID56031_Stress_FLASH_BEAT_epi_DS_AIF',  'NX20',   'GTPrep_2DT_Perf_AIF_2E_Lin_Mapping_MBF_MBV_Mask_AI_CMR_View_OFFLINE.xml', 'res',       'ref',   'IsmrmrdParameterMap_Siemens_Perfusion_NX20A.xsl'   ; ...            
    'SOLA/20201001_Marecav_pat1',           'meas_MID00577_FID56037_Rest_FLASH_BEAT_epi_DS_AIF',  'NX20',   'GTPrep_2DT_Perf_AIF_2E_Lin_Mapping_MBF_MBV_Mask_AI_CMR_View_OFFLINE.xml', 'res',       'ref',   'IsmrmrdParameterMap_Siemens_Perfusion_NX20A.xsl'   ; ...                
    
    'SOLA/loc_XA20',           'loc_XA20',  'NX20',   'Generic_Cartesian_Grappa.xml', 'res',       'ref',   'IsmrmrdParameterMap_Siemens_NX20A.xsl'   ; ...
    
    'perfusion/20200229_working_case',           'meas_MID00159_FID39583_SSFP_BEAT_epi_DS_AIF',  'NX20',   'GTPrep_2DT_Perf_AIF_2E_Lin_Mapping_MBF_MBV_Mask_AI_CMR_View_OFFLINE.xml', 'res',       'ref',   'IsmrmrdParameterMap_Siemens_Perfusion_NX20A.xsl'   ; ...            
    'perfusion/20200229_failed_case',           'meas_MID00488_FID39912_SSFP_BEAT_epi_DS_AIF',  'NX20',   'GTPrep_2DT_Perf_AIF_2E_Lin_Mapping_MBF_MBV_Mask_AI_CMR_View_OFFLINE.xml', 'res',       'ref',   'IsmrmrdParameterMap_Siemens_Perfusion_NX20A.xsl'   ; ...            
};

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end

if(nargout<=0)
    UTCases = [];
end


