function [perf, ori, moco, moco_norm, PD, input_for_filter, filtered, aif_acq_time, perf_acq_time, dst_acq_time, ... 
    aif_im, aif_moco, aif_cin, aif_cin_Gd, aif_cin_Gd_without_R2Star, aif_cin_Gd_baseline_corrected, ... 
    aif_cin_all_echo0_signal, aif_cin_all_echo1_signal, aif_cin_all_echo0_signal_after_R2StarCorrection, ...
    aif_cin_all_echo0_OverPD_after_R2StarCorrection, aif_cin_all_R2Star,  aif_cin_all_R2Star_SLEP, ... 
    aif_PD, aif_mask, aif_mask_final, aif_LV_mask_plot, aif, aif_baseline_corrected, ... 
    flow, Ki, PS, Vp, Visf, E, SDMap, Delay, ...
    BTEX_Flow_all, BTEX_PS_all, BTEX_Vp_all, BTEX_Visf_all, BTEX_cost_all, BTEX_flow_SD_all, Fermi_Delay] = read_in_GT_Perf_DebugOutput_results(resDir)
% read in Gadgetron perfusion debug output results

    slc = 0;   
    for n=1:8
        
        filename = ['perf_moco_upsampled_' num2str(n-1) '.hdr'];
        if(~isFileExist(fullfile(resDir, 'DebugOutput', filename)))
            break;
        end        
        slc = slc + 1;
    end
    
    disp(['Total ' num2str(slc) ' is found ...']);
    
    try
%         perf = load_array(resDir, 'CASignal_Perf_', slc);        
        perf = load_array2(resDir, 'PerfFlowMapping_Job_', slc, '_perf_moco_upsampled');        
        perf = permute(perf, [1 2 4 3]);
        perf = flipdim(perf, 2);
    catch
        perf = readGTPlusExportImageSeries_Squeeze(resDir, 111);
    end

    aif_acq_time = analyze75read(fullfile(resDir, 'DebugOutput', 'AIF_AcqTimes_0'));   
    perf_acq_time = load_array(resDir, 'Perf_AcqTimes_', slc);
    dst_acq_time = analyze75read(fullfile(resDir, 'DebugOutput', 'dstAcqTimes_0'));   
    
    ori = load_array(resDir, 'perf_', slc);   
    ori = flipdim(ori, 2);

    moco = load_array(resDir, 'input_for_SRNorm_', slc);   
    moco = flipdim(moco, 2);
    
    moco_norm = load_array(resDir, 'SRNorm_', slc);   
    moco_norm = flipdim(moco_norm, 2);

    PD = load_array(resDir, 'PD_', slc);   
    PD = flipdim(PD, 2);

    try
        input_for_filter = load_array(resDir, 'input_spatiotemporal_filter__row', slc);   
        input_for_filter = flipdim(input_for_filter, 2);
    catch
        input_for_filter = [];
    end
    
    try
        filtered = load_array(resDir, 'output_spatiotemporal_filter__row', slc);   
        filtered = flipdim(filtered, 2);
    catch
        filtered = [];
    end
    
    try
        r1 = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_moco.hdr'));
        r2 = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_moco_second_echo.hdr'));
        aif_moco = flipdim(cat(4, r1, r2), 2);
    catch
        aif_moco = [];
    end

    try
%         aif_im = readGTPlusExportImageSeries_Squeeze(resDir, 1104);
        aif_im = aif_moco;
    catch
        aif_im = [];
    end
    
    try
        aif_cin = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin.hdr'));
        aif_cin_Gd = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo0_LUTCorrection.hdr'));
        aif_cin_Gd_without_R2Star = analyze75read(fullfile(resDir, 'DebugOutput', 'cin_all_echo0_without_R2Star_LUTCorrection.hdr'));
        aif_cin_Gd_baseline_corrected = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_echo0_all_signal_baseline_corrected.hdr'));
        aif_cin_all_echo0_signal = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo0_signal.hdr'));
        aif_cin_all_echo1_signal = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo1_signal.hdr'));
        aif_cin_all_echo0_signal_after_R2StarCorrection = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo0_signal_after_R2StarCorrection.hdr'));
        aif_cin_all_echo0_OverPD_after_R2StarCorrection = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo0_OverPD_after_R2StarCorrection.hdr'));
        aif_cin_all_R2Star = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_R2Star.hdr'));
        aif_cin_all_R2Star_SLEP = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_R2Star_SLEP.hdr'));
        aif_PD = analyze75read(fullfile(resDir, 'DebugOutput', 'aifPD_for_TwoEcho_T2StartCorrection_0.hdr'));
        aif_mask = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_LV_mask_for_TwoEcho_T2StartCorrection_0.hdr'));
        aif_mask_final = analyze75read(fullfile(resDir, 'DebugOutput', 'AifLVMask_after_Picking.hdr'));

        aif_mask = flipdim(aif_mask, 2);
        aif_mask_final = flipdim(aif_mask_final, 2);
        aif_PD = flipdim(aif_PD, 2);

        aif_LV_mask_plot = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_LV_mask_plot_.hdr'));
        
        aif = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo0_LUTCorrection.hdr'));
    catch
        aif = [];
    end

    try
        aif_baseline_corrected = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_echo0_all_signal_baseline_corrected.hdr'));
    catch
        aif_baseline_corrected = [];
    end

    try
        flow = load_array(resDir, 'flow_maps_after_hole_filling_', slc);
        flow = squeeze(flow(:,:,end,:));
        flow = flipdim(flow, 2);
    catch
        flow = [];
    end

    try
        Ki = load_array(resDir, 'Ki_maps_after_hole_filling_', slc);        
        Ki = flipdim(Ki, 2);
    catch
        Ki = [];
    end

    try
        PS = load_array(resDir, 'PS_maps_after_hole_filling_', slc);
        PS = squeeze(PS(:,:,end,:));
        PS = flipdim(PS, 2);
    catch
        PS = [];
    end

    try
        Vp = load_array(resDir, 'blood_volume_maps_after_hole_filling_', slc);
        Vp = squeeze(Vp(:,:,end,:));
        Vp = flipdim(Vp, 2);
    catch
        Vp = [];
    end

    try
        Visf = load_array(resDir, 'interstitial_volume_maps_', slc);
        Visf = squeeze(Visf(:,:,end,:));
        Visf = flipdim(Visf, 2);
    catch
        Visf = [];
    end

    try
        E = load_array(resDir, 'E_maps_', slc);
        E = squeeze(E(:,:,end,:));
        E = flipdim(E, 2);
    catch
        E = [];
    end

    try
        for n=1:slc
            filename = ['BTEX_SD_maps_' num2str(n-1) '_0.hdr'];
            SDMap(:,:,n) = analyze75read(fullfile(resDir, 'DebugOutput', filename));
        end
        SDMap = flipdim(SDMap, 2);
    catch
        SDMap = [];
    end

    try
        Delay = load_array(resDir, 'BTEX_res_', slc);
        Delay = squeeze(Delay(:,:,end,:));
        Delay = flipdim(Delay, 2);
    catch
        Delay = [];
    end

    try
        Fermi_Delay = load_array(resDir, 'Fermi_res_', slc);
        Fermi_Delay = squeeze(Fermi_Delay(:,:,end,:));
        Fermi_Delay = flipdim(Fermi_Delay, 2);
    catch
        Fermi_Delay = [];
    end
    
    try
        BTEX_Flow_all = load_array(resDir, 'BTEX_Flow_all_', slc);
        BTEX_Flow_all = flipdim(BTEX_Flow_all, 2);
    catch
        BTEX_Flow_all = [];
    end

    try
        BTEX_PS_all = load_array(resDir, 'BTEX_PS_all_', slc);
        BTEX_PS_all = flipdim(BTEX_PS_all, 2);
    catch
        BTEX_PS_all = [];
    end

    try
        BTEX_Visf_all = load_array(resDir, 'BTEX_Visf_all_', slc);
        BTEX_Visf_all = flipdim(BTEX_Visf_all, 2);
    catch
        BTEX_Visf_all = [];
    end

    try
        BTEX_Vp_all = load_array(resDir, 'BTEX_Vp_all_', slc);
        BTEX_Vp_all = flipdim(BTEX_Vp_all, 2);
    catch
        BTEX_Vp_all = [];
    end

    try
        BTEX_cost_all = load_array(resDir, 'BTEX_cost_all_', slc);
        BTEX_cost_all = flipdim(BTEX_cost_all, 2);
    catch
        BTEX_cost_all = [];
    end

    try
        BTEX_flow_SD_all = load_array(resDir, 'BTEX_flow_SD_all_', slc);
        BTEX_flow_SD_all = flipdim(BTEX_flow_SD_all, 2);
    catch
        BTEX_flow_SD_all = [];
    end
end

function v = load_array(resDir, name, slc)
    try
        for n=1:slc
            filename = [name num2str(n-1) '.hdr'];
            v(:,:,:,n) = analyze75read(fullfile(resDir, 'DebugOutput', filename));
        end
    catch
        for n=1:slc
            filename = [name num2str(n-1) '_MAG.hdr'];
            v(:,:,:,n) = analyze75read(fullfile(resDir, 'DebugOutput', filename));
        end
    end
end

function v = load_array2(resDir, name, slc, name_after_slc)
    try
        for n=1:slc
            filename = [name num2str(n-1) name_after_slc '.hdr'];
            v(:,:,:,n) = analyze75read(fullfile(resDir, 'DebugOutput', filename));
        end
    catch
        for n=1:slc
            filename = [name num2str(n-1) name_after_slc '_MAG.hdr'];
            v(:,:,:,n) = analyze75read(fullfile(resDir, 'DebugOutput', filename));
        end
    end
end
