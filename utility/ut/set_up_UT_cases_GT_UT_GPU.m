function UTCases = set_up_UT_cases_GT_UT_GPU
% run the gt recon

% sub folder, data name, VB or VD, config xml file, result_folder, ground-truth folder
UTCases = {
'GT_UT/GPU/rtgrappa',  'acc_data_with_device_2', 'VB',   'grappa_float.xml',      'res',           'ref',   'IsmrmrdParameterMap.xsl'   ; ...   
'GT_UT/GPU/rtgrappa',  'acc_data_with_device_2', 'VB',   'grappa_float_cpu.xml',      'res',           'ref',   'IsmrmrdParameterMap.xsl'   ; ...   
'GT_UT/GPU/radial_phantom',  'meas_MID00133_FID20080_CV_Radial_Fixed_Angle_128_x8_32phs', 'VB',   'fixed_radial_mode1_realtime.xml',      'res',           'ref',   'IsmrmrdParameterMap.xsl' 
'GT_UT/snr_unit_recon_tpat3',           'meas_MID00171_FID02908_GRE_reps=150_TPAT=3',                   'VD',   'grappa_float.xml', 'res_gpu',       'ref_gpu',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...        
'GT_UT/snr_unit_recon_tpat3',           'meas_MID00171_FID02908_GRE_reps=150_TPAT=3',                   'VD',   'grappa_float_cpu.xml', 'res_cpu',       'ref_cpu',   'IsmrmrdParameterMap_Siemens.xsl'   ; ...        
};

N = size(UTCases, 1);

for i=1:N
    disp([num2str(i) '  -  ' UTCases{i, 1} '         -          ' UTCases{i, 2} ' - ' UTCases{i, 3} ' - ' UTCases{i, 4}]);
end
