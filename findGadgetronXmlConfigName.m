function configName = findGadgetronXmlConfigName(dataDir, datName)
% find the gadgetron xml config name
% first, it will try to find  the xml_raw.xml for tICEProgramName
% if this does not exist, it will try to determine it from config_buffer.xprot

configName = [];
        
xmlName = fullfile(dataDir, 'xml_raw.xml');

tICEProgramName = [];
if(isFileExist(xmlName))
    
    try
        ss = xml_load(fullfile(dataDir, 'xml_raw.xml'));       
        tICEProgramName = ss.YAPS(1).tICEProgramName  
        
        ind = strfind(tICEProgramName, 'IceProgramGadgetron_');
        if(~isempty(ind))
            tICEProgramName = tICEProgramName(ind(1):end)
        end
    catch
        tICEProgramName = [];
    end
end

if(isempty(tICEProgramName))
    % find ice prog name from config_buffer.xprot
    protName = fullfile(dataDir, 'config_buffer.xprot');

    fid = fopen(protName, 'r');
    prot = fscanf(fid,'%s');
    fclose(fid);
    
    ind = strfind(prot, '<ParamString."tICEProgramName">');
    if(~isempty(ind))
        prot2 = prot(ind(1):end);
        ind2 = strfind(prot2, 'IceProgramGadgetron_');
        if(~isempty(ind2))
            prot_end = prot2(ind2(1):end);
            ind_end = strfind(prot_end, '"');
            tICEProgramName = prot_end(1:ind_end(1)-1)
        end
    end
end

if(~isempty(strfind(tICEProgramName, 'IceProgramGadgetron')))
    switch tICEProgramName
        case {'IceProgramGadgetron_2DT'}
            configName = 'Generic_Cartesian_Grappa.xml';

        case {'IceProgramGadgetron_2DT_GFactor'}
            configName = 'Generic_Cartesian_Grappa_SNR.xml';

        case {'IceProgramGadgetron_2DT_RealTimeCine'}
            configName = 'GTPrep_2DT_RealTimeCine_PF_Handling_PhysioInterp.xml';

        case {'IceProgramGadgetron_2DT_RealTimeCine_KSpaceBinning'}
            configName = 'GTPrep_2DT_RTCine_KspaceBinning.xml';

        case {'IceProgramGadgetron_2DT_RealTimeCine_KSpaceBinning_Amazon_Cloud'}
            configName = 'GTPrep_2DT_RTCine_KspaceBinning_SingleLayer_Gateway.xml';

        case {'IceProgramGadgetron_2DT_RealTimeCine_NonLinear'}
            configName = 'GTPrep_2DT_RTCine_L1SPIRIT_SLEP_PhysioInterp.xml';

        case {'IceProgramGadgetron_2DT_RealTimeCine_NonLinear_Amazon_Cloud'}
            configName = 'GTPrep_2DT_RTCine_L1SPIRIT_SLEP_PhysioInterp_DualLayer_Gateway.xml';

        case {'IceProgramGadgetron_2DT_RetroGated_Cine'}
            configName = 'GTPrep_2DT_RetroGated_Cine.xml';

        case {'IceProgramGadgetron_2DT_RetroGated_Flow'}
            configName = 'GTPrep_2DT_RetroGated_Flow.xml';

        case {'IceProgramGadgetron_2DT_RetroGated_Cine_L1SPIRIT_SLEP'}
            configName = 'GTPrep_2DT_RetroGated_Cine_L1SPIRIT_SLEP.xml';

        case {'IceProgramGadgetron_Perfusion_Linear_AIF_2echo_TPAT2'}
            configName = 'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_QuantitativeFlow_Mapping.xml';

        case {'IceProgramGadgetron_Perfusion_Linear_AIF_2echo_TPAT2_NonLinear'}
            configName = 'GTPrep_2DT_Perfusion_AIF_TwoEchoes_Interleaved_R2_QuantitativeFlow_Mapping_NonLinear.xml';

        case {'IceProgramGadgetron_PSIR'}
            configName = 'GTPrep_2DT_LGE_MOCO_AVE_OnTheFly.xml';

        case {'IceProgramGadgetron_DB_PSIR'}
            configName = 'GTPrep_2DT_LGE_MOCO_AVE_OnTheFly.xml';

        case {'IceProgramGadgetron_T2W'}
            configName = 'GTPrep_2DT_T2W_MOCO_AVE_OnTheFly.xml';

        case {'IceProgramGadgetron_T2StarMapping'}
            configName = 'GTPrep_2DT_T2Star_Mapping.xml';                                   

        otherwise
            configName = [];
    end
else
    nameUsed = lower(datName);
    if(~isempty(strfind(nameUsed, 'loc')) || ~isempty(strfind(nameUsed, 'localizer')) )
        configName = 'Generic_Cartesian_Grappa.xml';
    elseif (~isempty(strfind(nameUsed, 'realtime')) || ~isempty(strfind(nameUsed, 'realtime_gt')) )
        configName = 'GTPrep_2DT_RealTimeCine_PF_Handling_PhysioInterp.xml';
    elseif (~isempty(strfind(nameUsed, 'fisp_retro')) )
        configName = 'GTPrep_2DT_RetroGated_Cine.xml';
    elseif (~isempty(strfind(nameUsed, 'molli')) )
        configName = 'GTPrep_2DT_MOLLI.xml';
    elseif (~isempty(strfind(nameUsed, 'sasha')) )
        configName = 'GTPrep_2DT_SASHA.xml';
    elseif (~isempty(strfind(nameUsed, 'fw')) || ~isempty(strfind(nameUsed, 'fat'))  )
        configName = 'Generic_Cartesian_Grappa_FatWater.xml';
    elseif (~isempty(strfind(nameUsed, 't2_')) )
        configName = 'Generic_Cartesian_Grappa.xml';
        % configName = 'GTPrep_2DT_MOCO_AVE_T2_Mapping.xml';
    elseif (~isempty(strfind(nameUsed, 'tse')) )
        configName = 'Generic_Cartesian_Grappa.xml';
    elseif (~isempty(strfind(nameUsed, 'tag')) )
        configName = 'Generic_Cartesian_Grappa.xml';
    end
    
    if(isempty(configName))
        switch tICEProgramName
            case {'IceProgram_2DT'}
                configName = 'Generic_Cartesian_Grappa.xml';
            otherwise
        end
    end
end

if(isempty(configName))
    tProtocolName = lower(ss.MEAS(1).tProtocolName)
    
    if(strfind(tProtocolName, 'realtime_cine'))
        configName = 'GTPrep_2DT_RealTimeCine_PF_Handling_PhysioInterp.xml';
    end
    
    if(strfind(tProtocolName, 'retro'))
        configName = 'GTPrep_2DT_RetroGated_Cine.xml';
    end
end
        