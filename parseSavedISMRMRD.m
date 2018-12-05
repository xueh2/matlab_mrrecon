function [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name)
% [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD('LGE_MOCO_AVE_OnTheFly_41672_7432022_7432027_319_20160329-185042')

len = length(name);
ind = strfind(name, '_');
ind2 = strfind(name, '-');

if(isempty(ind))
    error(['Empty name ' name])
end

st = name(ind(end)+1:end);
if(isempty(strfind(st, ':')))
    study_time = name(ind2(end)+1:end);
    study_dates = name(ind(end)+1:ind2(end)-1);
else
    study_time = [name(len-7:len-6) name(len-4:len-3) name(len-1:len)];
    study_dates = [name(len-18:len-15) name(len-13:len-12) name(len-9:len-10)];
end

study_year = study_dates(1:4);
study_month = study_dates(5:6);
study_day = study_dates(7:8);

name2 = name(1:ind(end)-1);
ind = strfind(name2, '_');

measurementID = name2(ind(end)+1:end);
studyID = name2(ind(end-1)+1:ind(end)-1);
patientID = name2(ind(end-2)+1:ind(end-1)-1);
scannerID = name2(ind(end-3)+1:ind(end-2)-1);

xmlString = name2(1:ind(end-3)-1);

configName = [];
if(~isempty(strfind(xmlString, 'ISMRMRD_Noise_dependency')))
    configName = 'default_measurement_dependencies.xml';
    
elseif(~isempty(strfind(xmlString, 'DB_LGE_MOCO_AVE_OnTheFly')))
    configName = 'GTPrep_2DT_DB_LGE_MOCO_AVE_OnTheFly_dstore.xml';
    %configName = 'GTPrep_2DT_DB_LGE_MOCO_AVE_OnTheFly_MoreOutputs_dstore.xml';
    
elseif(~isempty(strfind(xmlString, 'LGE_MOCO_AVE')))
    configName = 'GTPrep_2DT_LGE_MOCO_AVE_OnTheFly_dstore.xml';
    %configName = 'GTPrep_2DT_LGE_MOCO_AVE_OnTheFly_MoreOutputs_dstore.xml';
    
elseif(~isempty(strfind(xmlString, 'T2W')))
    configName = 'GTPrep_2DT_T2W_MOCO_AVE_OnTheFly_istore.xml';
    
elseif(~isempty(strfind(xmlString, 'T2Star_Mapping')))
    configName = 'GTPrep_2DT_T2Star_Mapping_dstore.xml';    
    
elseif(~isempty(strfind(xmlString, 'T2_Mapping')))
    configName = 'GTPrep_2DT_MOCO_AVE_T2_Mapping_dstore.xml';
    
elseif(~isempty(strfind(xmlString, 'Perfusion_AIF_TwoEchoes_Interleaved')))
    configName = 'GTPrep_2DT_Perf_AIF_2E_Lin_Mapping_OFFLINE_dstore.xml';
elseif(~isempty(strfind(xmlString, 'Perfusion_AIF_2E_Lin_Cloud')))
    configName = 'GTPrep_2DT_Perf_AIF_2E_Lin_Mapping_OFFLINE_dstore.xml';
elseif(~isempty(strfind(xmlString, 'Perfusion_AIF_2E_NL_Cloud')))
    configName = 'GTPrep_2DT_Perf_AIF_2E_NL_Mapping_MBF_MBV_Mask_Gateway.xml';
 elseif(~isempty(strfind(xmlString, 'Perfusion_AIFR3_2E')))
    configName = 'GTPrep_2DT_Perf_AIFR3_2E_Lin_Mapping_MBF_MBV_Mask_OFFLINE.xml';
   
elseif(~isempty(strfind(xmlString, 'Retro_Flow')))
    configName = 'GTPrep_2DT_RetroGated_Flow.xml';    
elseif(~isempty(strfind(xmlString, 'Retro_NLin_Flow')))
    configName = 'GTPrep_2DT_RetroGated_Flow_SLEP.xml'; 
    
elseif(~isempty(strfind(xmlString, 'Prospective_Cine_3D')))
    configName = 'GTPrep_3DT_RetroGated_Cine_istore.xml'; 
    
elseif(~isempty(strfind(xmlString, 'Cmr_Binning_Cloud')))
    configName = 'CMR_2DT_RTCine_KspaceBinning_Cloud.xml';         
elseif(~isempty(strfind(xmlString, 'Binning')))
    configName = 'CMR_2DT_RTCine_KspaceBinning.xml';     
    
elseif(~isempty(strfind(xmlString, 'Cmr_Cine_NL_Cloud')))
    configName = 'Generic_Cartesian_NonLinear_Spirit_RealTimeCine_Cloud.xml';         
elseif(~isempty(strfind(xmlString, 'Cine_NL')))
    configName = 'Generic_Cartesian_NonLinear_Spirit_RealTimeCine.xml';        
 
elseif(~isempty(strfind(xmlString, 'Retro_NLin_Cine')))
    configName = 'GTPrep_2DT_RetroGated_Cine_SLEP_Gateway.xml';        
elseif(~isempty(strfind(xmlString, 'Retro_Lin_Cine')))
    configName = 'Generic_2DT_RetroGated_Cine_ECG_dstore.xml';  
    
elseif(~isempty(strfind(xmlString, 'RT_Cine_LIN')))
    configName = 'Generic_RTCine_PInterp_Fil_ECG_dstore.xml';            
    
elseif(~isempty(strfind(xmlString, 'FatWater_Nav3D')))    
    configName = 'GTPrep_3DT_FatWater_Diego.xml'; 
elseif(~isempty(strfind(xmlString, 'FatWaterCine')))
    configName = 'GTPrep_2DT_FW_CINE_Diego.xml'; 
elseif(~isempty(strfind(xmlString, 'FatWater_MOCO_AVE_PSIR')))
    configName = 'GTPrep_2DT_FW_MOCO_AVE_PSIR_Diego_AVE.xml'; 
elseif(~isempty(strfind(xmlString, 'FatWater_MOCO_AVE')))
    configName = 'GTPrep_2DT_FW_MOCO_AVE_Diego.xml';     
elseif(~isempty(strfind(xmlString, 'FatWater')))
    configName = 'GTPrep_2DT_FW_Diego.xml';         

end