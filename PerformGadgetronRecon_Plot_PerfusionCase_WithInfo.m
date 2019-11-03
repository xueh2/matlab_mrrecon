
function h_flow = PerformGadgetronRecon_Plot_PerfusionCase_WithInfo(resDir, restDir, scannerID, patientID, studyID, study_dates, study_time, onlyReview, checkProcessed, baseDir)
% h_flow = PerformGadgetronRecon_Plot_PerfusionCase_WithInfo(restDir, scannerID, patientID, studyID, study_dates, onlyReview, baseDir)

if(nargin < 7)
    onlyReview = 0;
end

if(nargin < 8)
    checkProcessed = 1;
end

if(nargin < 9)
    baseDir = resDir;
end

h_flow = 0;

resDir = [scannerID '_' patientID '_' studyID '_' study_dates '_' study_time];

figDir = fullfile(baseDir, study_dates, ['Perfusion_AIFR3_2E_Interleaved_' resDir '_Figure']);

if(checkProcessed & isFileExist(fullfile(figDir, 'rest.mat')))
    disp(['Already processed - ' restDir]);
    h_flow = -1;
    return;
end

if(~onlyReview)   
    try
        [perf, ori, moco, moco_norm, PD, input_for_filter, filtered, aif_acq_time, perf_acq_time, dst_acq_time, ... 
    aif_im, aif_moco, aif_cin, aif_cin_Gd, aif_cin_Gd_without_R2Star, aif_cin_Gd_baseline_corrected, ... 
    aif_cin_all_echo0_signal, aif_cin_all_echo1_signal, aif_cin_all_echo0_signal_after_R2StarCorrection, ...
    aif_cin_all_echo0_OverPD_after_R2StarCorrection, aif_cin_all_R2Star,  aif_cin_all_R2Star_SLEP, ... 
    aif_PD, aif_mask, aif_mask_final, aif_LV_mask_plot, aif, aif_baseline_corrected, aif_plots, ... 
    flow, Ki, PS, Vp, Visf, E, SDMap, Delay, ...
    BTEX_Flow_all, BTEX_PS_all, BTEX_Vp_all, BTEX_Visf_all, BTEX_cost_all, BTEX_flow_SD_all, BTEX_Tc_all, Fermi_Delay] = read_in_GT_Perf_DebugOutput_results(restDir);
                 
        has_rest = 1;
       
    catch
        has_rest = 0;
        perf = 0;
        fa = 0;
        a = 0;
        aif = 0;
    end  
end

% scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time

if(~onlyReview)
    if(has_rest)        
        if(exist(figDir)==7)
            rmdir(figDir, 's');
        end
        mkdir(figDir);
        
        if(has_rest)
            % save(fullfile(figDir, 'rest.mat'), 'perf', 'ori', 'moco', 'moco_norm', 'aif_im', 'aif_moco', 'aif', 'aif_baseline_corrected', 'aif_cin', 'aif_cin_Gd', 'aif_cin_Gd_without_R2Star', 'aif_cin_Gd_baseline_corrected', 'aif_cin_all_echo0_signal', 'aif_cin_all_echo1_signal', 'aif_cin_all_echo0_signal_after_R2StarCorrection', 'aif_cin_all_echo0_OverPD_after_R2StarCorrection', 'aif_cin_all_R2Star', 'aif_cin_all_R2Star_SLEP', 'aif_PD', 'aif_mask', 'aif_mask_final', 'flow', 'Ki', 'PS', 'Vp', 'Visf', 'E', 'SDMap', 'Delay');
            save(fullfile(figDir, 'rest.mat'), 'perf', 'ori', 'moco', 'moco_norm', 'PD', 'input_for_filter', 'filtered', 'aif_acq_time', 'perf_acq_time', 'dst_acq_time', ... 
                'aif_im', 'aif_moco', 'aif', 'aif_baseline_corrected', 'aif_cin', ... 
                'aif_cin_Gd', 'aif_cin_Gd_without_R2Star', 'aif_cin_Gd_baseline_corrected', ... 
                'aif_cin_all_echo0_signal', 'aif_cin_all_echo1_signal', 'aif_cin_all_echo0_signal_after_R2StarCorrection', ...
                'aif_cin_all_echo0_OverPD_after_R2StarCorrection', 'aif_cin_all_R2Star', 'aif_cin_all_R2Star_SLEP', ...
                'aif_PD', 'aif_mask', 'aif_mask_final', ...
                'flow', 'Ki', 'PS', 'Vp', 'Visf', 'E', 'SDMap', 'Delay', ...
                'BTEX_Flow_all', 'BTEX_PS_all', 'BTEX_Visf_all', 'BTEX_Vp_all', 'BTEX_cost_all', 'BTEX_flow_SD_all', 'BTEX_Tc_all');
        end
    end
else
    has_rest = 1;
end

if(has_rest)
    if(~onlyReview)
        slc = size(Ki, 4);
        m = size(Ki, 3);

        scalingFactor = 10;

        figure; imagescn(Ki, [0 6], [m slc], scalingFactor); PerfColorMap;
    end
    scrsz = get(0, 'ScreenSize');
    
    figName = fullfile(figDir, [resDir '_Rest_AIF_Acq_Plot.fig']);
    if(onlyReview)
        if(isFileExist([figName '.fig']))
            openfig(figName);
        end
    else
        h = figure('visible', 'off', 'Name','AIF Acqusition Time Plot','NumberTitle','off', 'Position', [100 100 3*1024 1.5*768]);            
        hold on

        aif_v = interp1(dst_acq_time, aif_cin_Gd_baseline_corrected, aif_acq_time, 'linear');
        aif_t = (aif_acq_time-dst_acq_time(1))/1000;
        plot(aif_t, aif_v, 'k+', 'MarkerSize',12, 'LineWidth', 2);

        perf_acq_time = squeeze(perf_acq_time);

        C = {'m','g','r', 'k', 'c', 'b', 'y', 'b'};
        for ss=1:size(perf_acq_time, 2)
            perf_t = interp1(dst_acq_time, aif_cin_Gd_baseline_corrected, perf_acq_time(:,ss), 'linear');            
            plot((perf_acq_time(:,ss)-dst_acq_time(1))/1000, perf_t, [C{ss} '.'], 'MarkerSize',16);
        end
        xlabel('Acqusition time, in s');
        ylabel('AIF Gd, mmol/L')
        legend('AIF', 'Perf slice 1', 'Perf slice 2', 'Perf slice 3')
        title('Stress')

        xlim_v = [min(aif_t)-0.5 max(aif_t)+0.5];
        ylim_v = [min(aif_cin_Gd_baseline_corrected)-1 max(aif_cin_Gd_baseline_corrected)+1];

        plot((dst_acq_time-dst_acq_time(1))/1000, aif_cin_Gd_baseline_corrected, 'b--');            

        plot([xlim_v(1) xlim_v(2)], [ylim_v(1)+0.5 ylim_v(1)+0.5], 'b-')
        plot(aif_t, ylim_v(1)+0.5, 'k+', 'MarkerSize',12, 'LineWidth', 2)

        xlim(xlim_v);
        ylim(ylim_v);
        hold off
        box on
        grid on
        saveas(h, figName, 'fig')
    end
    
    % ----------------------
    
    figName = fullfile(figDir, [resDir '_Rest_Ki_TwoComp.fig']);
    if(onlyReview)
        if(isFileExist([figName '.fig']))
            openfig(figName);
        end
    else
        if(m==4)
            h = figure('visible', 'off', 'Name','Ki Maps - TwoCompExp','NumberTitle','off'); imagescn(Ki(:,:,3,:), [0 6], [1 slc], scalingFactor); PerfColorMap;
            saveas(h, figName, 'fig')
        end
    end
    
    figName = fullfile(figDir, [resDir '_Rest_Gd.fig']);
    if(onlyReview)
        if(isFileExist([figName '.fig']))
            openfig(figName);
        end
    else
        h = figure('visible', 'off', 'Name','Rest Gd','NumberTitle','off'); imagescn(perf, [0 1.5], [1 slc], 10, 4);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_AIF_ORI' '.fig']);
    if(onlyReview)
        if(isFileExist([figName '.fig']))
            openfig(figName);
        end
    else        
        h = figure('visible', 'off', 'Name','AIF Original','NumberTitle','off'); imagescn(aif_im, [], [], 10, 4);
        saveas(h, figName, 'fig');
    end
            
    figName = fullfile(figDir, [resDir '_Rest_PDE_FlowMap.fig']);
    if(onlyReview)
        if(isFileExist([figName '.fig']))
            h_flow = openfig(figName);
        end
    else
        h = figure('visible', 'off', 'Name','Flow maps','NumberTitle','off'); imagescn(flow(:,:,:,end), [0 6], [1 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_PDE_Visf.fig']);
    if(onlyReview)
        if(isFileExist([figName '.fig']))
            openfig(figName);
        end
    else
        if(~isempty(Visf))
            h = figure('visible', 'off', 'Name','PDE Visf','NumberTitle','off'); imagescn(Visf(:,:,:,end), [0 80], [1 slc], scalingFactor); ECVColorMap;
            saveas(h, figName, 'fig')
        end
    end
    
    figName = fullfile(figDir, [resDir '_Rest_PDE_PS.fig']);
    if(onlyReview)
        if(isFileExist([figName '.fig']))
            openfig(figName);
        end
    else
        if(~isempty(PS))
            h = figure('visible', 'off', 'Name','PDE PS','NumberTitle','off'); imagescn(PS(:,:,:,end), [0 2], [1 slc], scalingFactor); PSColorMap;
            saveas(h, figName, 'fig')
        end
    end
    
    figName = fullfile(figDir, [resDir '_Rest_PDE_E.fig']);
    if(onlyReview)
        if(isFileExist([figName '.fig']))
            openfig(figName);
        end
    else
        if(~isempty(E))
            h = figure('visible', 'off', 'Name','PDE E','NumberTitle','off');; imagescn(E(:,:,:,end), [0 2], [1 slc], scalingFactor); PerfColorMap;
            saveas(h, figName, 'fig')
        end
    end
    
    figName = fullfile(figDir, [resDir '_Rest_PDE_Vp.fig']);
    if(onlyReview)
        if(isFileExist([figName '.fig']))
            openfig(figName);
        end
    else
        if(~isempty(Vp))
            h = figure('visible', 'off', 'Name','PDE Vb','NumberTitle','off');; imagescn(Vp(:,:,:,end), [0 20], [1 slc], scalingFactor); MBVColorMap;
            saveas(h, figName, 'fig')
        end
    end
    
    figName = fullfile(figDir, [resDir '_Rest_PDE_Tc.fig']);
    if(onlyReview)
        if(isFileExist([figName '.fig']))
            openfig(figName);
        end
    else
        if(~isempty(BTEX_Tc_all))
            h = figure('visible', 'off', 'Name','PDE Tc','NumberTitle','off'); imagescn(BTEX_Tc_all(:,:,:,end), [0 12], [1 slc], scalingFactor); MBVColorMap;
            saveas(h, figName, 'fig')
        end
    end
    
    try
        bulls_eye_rest = perf_load_Bulls_eye_debug_data(fullfile(restDir, 'DebugOutput'), 1, 0);
        [h_fmap_seg, h_Gd_seg, h_mask_seg] = perf_plot_loaded_Bulls_eye_debug_data(bulls_eye_rest);

        figName = fullfile(figDir, [resDir '_Rest_fmap_with_SectorContours' '.fig']);
        saveas(h_fmap_seg, figName, 'fig');
        figName = fullfile(figDir, [resDir '_Rest_Gd_with_SectorContours' '.fig']);
        saveas(h_Gd_seg, figName, 'fig');
        figName = fullfile(figDir, [resDir '_Rest_SectorMasks' '.fig']);
        saveas(h_mask_seg, figName, 'fig');
        save(fullfile(figDir, 'bulls_eye_rest'), 'bulls_eye_rest');
    catch
    end
            
    figName = fullfile(figDir, [resDir '_Rest_Ori.fig']);
    if(onlyReview)
        if(isFileExist([figName '.fig']))
            openfig(figName);
        end
    else
        h = figure('visible', 'off', 'Name','Rest Original','NumberTitle','off'); imagescn(ori, [], [1 slc], scalingFactor, 3);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_MOCO.fig']);
    if(onlyReview)
        if(isFileExist([figName '.fig']))
            openfig(figName);
        end
    else
        h = figure('visible', 'off', 'Name','Rest MOCO','NumberTitle','off'); imagescn(moco, [], [1 slc], scalingFactor, 3);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_MOCO_NORM.fig']);
    if(onlyReview)
        if(isFileExist([figName '.fig']))
            openfig(figName);
        end
    else
        h = figure('visible', 'off', 'Name','Rest MOCO NORM','NumberTitle','off'); imagescn(moco_norm, [], [1 slc], scalingFactor, 3);
        saveas(h, figName, 'fig')
    end
    
    if(~isempty(input_for_filter))
        figName = fullfile(figDir, [resDir '_Rest_MOCO_FIL.fig']);
        if(onlyReview)
            if(isFileExist([figName '.fig']))
                openfig(figName);
            end
        else
            h = figure('visible', 'off', 'Name','Rest MOCO Filtered','NumberTitle','off'); imagescn(cat(4, input_for_filter, filtered), [], [2 slc], scalingFactor, 3);
            saveas(h, figName, 'fig')
        end
    end
    
    NN = size(perf, 4);
    
    figName = fullfile(figDir, [resDir '_Rest_CoMOCO' '.fig']);
    if(onlyReview & isFileExist(figName))
        if(isFileExist(figName) ) 
            openfig(figName);
        else
            moco_pd = cat(3, PD, moco);
            rr = moco_pd(:,:,[1 4 round(NN/2) NN],:);

            h = figure('visible', 'off', 'Name',[ resDir '_Rest SR-PD CoMOCO'],'NumberTitle','off'); imagescn(rr, [], [slc 4], scalingFactor);
            saveas(h, figName, 'fig');
        end
    else
        moco_pd_rest = cat(3, PD, moco);
        rr = moco_pd_rest(:,:,[1 size(PD, 4)+1 end],:);

        h = figure('visible', 'off', 'Name',[ resDir '_Rest SR-PD CoMOCO'],'NumberTitle','off'); imagescn(rr, [], [slc 3], scalingFactor);
        saveas(h, figName, 'fig');
    end
            
    figName = fullfile(figDir, [resDir '_Rest_AIF_MOCO.fig']);
    if(onlyReview)
        if(isFileExist([figName '.fig']))
            openfig(figName);
        end
    else
        h = figure('visible', 'off', 'Name','AIF MOCO','NumberTitle','off'); imagescn(aif_moco, [], [], 10, 3);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_PDE_Delay' '.fig']);
    if(onlyReview)
        if(isFileExist([figName '.fig']))
            openfig(figName);
        end
    else
        if(~isempty(Delay))
            h = figure('visible', 'off', 'Name','PDE Delay','NumberTitle','off');; imagescn(Delay(:,:,:,end), [0 8], [1 slc], scalingFactor);
            saveas(h, figName, 'fig');
        end
    end
        
    delta = 0.5;

    figName = fullfile(figDir, [resDir '_AIF_Rest_Curves.fig']);

    if(onlyReview)
        if(isFileExist([figName '.fig']))
            openfig(figName);
        end
    else        
        h = figure
        hold on
        plot(delta*[0:numel(aif)-1], aif, 'LineWidth',2)
        hold off
        set(gcf, 'Units', 'normalized', 'Position', [0.2 0.2 0.5 0.5])
        box on
        grid on
        grid MINOR
        legend('rest')
        xlabel('second')
        ylabel('Gd [mmol/ml]')
        title('AIF in Gd')

        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_AIF_Rest_Curves' '.fig']);
    if(onlyReview)
        if(isFileExist([figName '.fig']))
            openfig(figName);
        end
    else        
        h = figure
        hold on
        plot(delta*[0:numel(aif_cin_all_echo0_signal)-1], aif_cin_all_echo0_signal, 'b', 'LineWidth',2);
        plot(delta*[0:numel(aif_cin_all_echo1_signal)-1], aif_cin_all_echo1_signal, 'k', 'LineWidth',2)
        plot(delta*[0:numel(aif_cin_all_echo0_signal_after_R2StarCorrection)-1], aif_cin_all_echo0_signal_after_R2StarCorrection, 'r', 'LineWidth',2)
        hold off
        set(gcf, 'Units', 'normalized', 'Position', [0.2 0.2 0.5 0.5])
        box on
        grid on
        grid MINOR
        legend('first echo', 'second echo', 'first echo with T2* correction')
        xlabel('second')
        ylabel('internsity')
        title('AIF T2* correction, rest')

        saveas(h, figName, 'fig');
    end 
    
    figName = fullfile(figDir, [resDir '_Rest_AIF_Mask' '.fig']);
    if(onlyReview)
        if(isFileExist([figName '.fig']))
            openfig(figName);
        end
    else        
                
        h = figure('visible', 'off', 'Name','AIF MOCO with Mask, Rest','NumberTitle','off'); imagescn(aif_LV_mask_plot, [], [], 10);
        saveas(h, figName, 'fig');
    end
end
