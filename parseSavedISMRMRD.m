function [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(name)
% [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD('LGE_MOCO_AVE_OnTheFly_41672_7432022_7432027_319_20160329-185042')

len = length(name);
ind = strfind(name, '_');

if(isempty(ind))
    error(['Empty name ' name])
end

st = name(ind(end)+1:end);
if(isempty(strfind(st, ':')))
    study_time = name(len-5:len);
    study_dates = name(len-14:len-7);
else
    study_time = [name(len-7:len-6) name(len-4:len-3) name(len-1:len)];
    study_dates = [name(len-18:len-15) name(len-13:len-12) name(len-9:len-10)];
end

study_year = study_dates(1:4);
study_month = study_dates(5:6);
study_day = study_dates(7:8);

name2 = name(1:len-16);
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
    
elseif(~isempty(strfind(xmlString, 'LGE_MOCO_AVE')))
    configName = 'GTPrep_2DT_LGE_MOCO_AVE_OnTheFly_dstore.xml';
    
elseif(~isempty(strfind(xmlString, 'T2W')))
    configName = 'GTPrep_2DT_T2W_MOCO_AVE_OnTheFly_istore.xml';
    
elseif(~isempty(strfind(xmlString, 'T2Star_Mapping')))
    configName = 'GTPrep_2DT_T2Star_Mapping_dstore.xml';    
    
elseif(~isempty(strfind(xmlString, 'Perfusion_AIF_TwoEchoes_Interleaved')))
    configName = 'GTPrep_2DT_Perf_AIF_2E_Lin_Mapping_OFFLINE_dstore.xml';
elseif(~isempty(strfind(xmlString, 'Perfusion_AIF_2E_Lin_Cloud')))
    configName = 'GTPrep_2DT_Perf_AIF_2E_Lin_Mapping_OFFLINE_dstore.xml';
elseif(~isempty(strfind(xmlString, 'Perfusion_AIF_2E_NL_Cloud')))
    configName = 'GTPrep_2DT_Perf_AIF_2E_Lin_Mapping_OFFLINE_dstore.xml';
    
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
    
elseif(~isempty(strfind(xmlString, 'FW')))
    configName = 'GTPrep_2DT_FatWater.xml';         

end