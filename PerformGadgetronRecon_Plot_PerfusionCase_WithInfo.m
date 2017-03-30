
function PerformGadgetronRecon_Plot_PerfusionCase_WithInfo(resDir, restDir, scannerID, patientID, studyID, study_dates, onlyReview, baseDir)
% PerformGadgetronRecon_Plot_PerfusionCase_WithInfo(restDir, scannerID, patientID, studyID, study_dates, onlyReview, baseDir)

if(nargin < 7)
    onlyReview = 0;
end

if(nargin < 8)
    baseDir = resDir;
end

if(~onlyReview)   
    try
        [rest_perf, ori_rest, moco_rest, moco_norm_rest, PD_rest, ... 
            aif_im_rest, aif_moco_rest, aif_rest_cin, aif_rest_cin_Gd, aif_rest_cin_Gd_without_R2Star, aif_rest_cin_Gd_baseline_corrected, ... 
            aif_rest_cin_all_echo0_signal, aif_rest_cin_all_echo1_signal, aif_rest_cin_all_echo0_signal_after_R2StarCorrection, ...
            aif_rest_cin_all_echo0_OverPD_after_R2StarCorrection, aif_rest_cin_all_R2Star,  aif_rest_cin_all_R2Star_SLEP, ... 
            aif_rest_PD, aif_rest_mask, aif_rest_mask_final, aif_rest, aif_rest_baseline_corrected, ... 
            flow_rest, Ki_rest, PS_rest, Vp_rest, Visf_rest, E_rest, SDMap_rest, Delay_rest, ...
            BTEX_Flow_all_rest, BTEX_PS_all_rest, BTEX_Vp_all_rest, BTEX_Visf_all_rest, BTEX_cost_all_rest, BTEX_flow_SD_all_rest] = read_in_GT_Perf_DebugOutput_results(restDir);
                 
        has_rest = 1;
       
    catch
        has_rest = 0;
        rest_perf = 0;
        fa = 0;
        a = 0;
        aif_rest = 0;
    end  
end

% scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time
resDir = [scannerID '_' patientID '_' studyID '_' study_dates];

figDir = fullfile(baseDir, study_dates, ['Perfusion_AIF_TwoEchoes_Interleaved_R2_' resDir '_Figure']);

if(~onlyReview)
    if(has_rest)        
        if(exist(figDir)==7)
            rmdir(figDir, 's');
        end
        mkdir(figDir);
        
        if(has_rest)
            % save(fullfile(figDir, 'rest.mat'), 'rest_perf', 'ori_rest', 'moco_rest', 'moco_norm_rest', 'aif_im_rest', 'aif_moco_rest', 'aif_rest', 'aif_rest_baseline_corrected', 'aif_rest_cin', 'aif_rest_cin_Gd', 'aif_rest_cin_Gd_without_R2Star', 'aif_rest_cin_Gd_baseline_corrected', 'aif_rest_cin_all_echo0_signal', 'aif_rest_cin_all_echo1_signal', 'aif_rest_cin_all_echo0_signal_after_R2StarCorrection', 'aif_rest_cin_all_echo0_OverPD_after_R2StarCorrection', 'aif_rest_cin_all_R2Star', 'aif_rest_cin_all_R2Star_SLEP', 'aif_rest_PD', 'aif_rest_mask', 'aif_rest_mask_final', 'flow_rest', 'Ki_rest', 'PS_rest', 'Vp_rest', 'Visf_rest', 'E_rest', 'SDMap_rest', 'Delay_rest');
            save(fullfile(figDir, 'rest.mat'), 'rest_perf', 'ori_rest', 'moco_rest', 'moco_norm_rest', ... 
                'aif_im_rest', 'aif_moco_rest', 'aif_rest', 'aif_rest_baseline_corrected', 'aif_rest_cin', ... 
                'aif_rest_cin_Gd', 'aif_rest_cin_Gd_without_R2Star', 'aif_rest_cin_Gd_baseline_corrected', ... 
                'aif_rest_cin_all_echo0_signal', 'aif_rest_cin_all_echo1_signal', 'aif_rest_cin_all_echo0_signal_after_R2StarCorrection', ...
                'aif_rest_cin_all_echo0_OverPD_after_R2StarCorrection', 'aif_rest_cin_all_R2Star', 'aif_rest_cin_all_R2Star_SLEP', ...
                'aif_rest_PD', 'aif_rest_mask', 'aif_rest_mask_final', ...
                'flow_rest', 'Ki_rest', 'PS_rest', 'Vp_rest', 'Visf_rest', 'E_rest', 'SDMap_rest', 'Delay_rest', ...
                'BTEX_Flow_all_rest', 'BTEX_PS_all_rest', 'BTEX_Visf_all_rest', 'BTEX_Vp_all_rest', 'BTEX_cost_all_rest', 'BTEX_flow_SD_all_rest');
        end
    end
else
    has_rest = 1;
end

if(has_rest)
    if(~onlyReview)
        slc = size(Ki_rest, 4);
        m = size(Ki_rest, 3);

        scalingFactor = 10;

        figure; imagescn(Ki_rest, [0 6], [m slc], scalingFactor); PerfColorMap;
    end
    scrsz = get(0, 'ScreenSize');
    
    figName = fullfile(figDir, [resDir '_Rest_Ki_TwoComp']);
    if(onlyReview)
        openfig(figName);
    else
        if(m==4)
            h = figure('Name','Ki Maps - TwoCompExp','NumberTitle','off'); imagescn(Ki_rest(:,:,3,:), [0 6], [1 slc], scalingFactor); PerfColorMap;
            saveas(h, figName, 'fig')
        end
    end
    
    figName = fullfile(figDir, [resDir '_Rest_Gd']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','Rest Gd','NumberTitle','off'); imagescn(rest_perf, [0 1.5], [1 slc], 10, 4);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_AIF_ORI' '.fig']);
    if(onlyReview)
        openfig(figName);
    else        
        h = figure('Name','AIF Original','NumberTitle','off'); imagescn(aif_im_rest, [], [], 10, 4);
        saveas(h, figName, 'fig');
    end
            
    figName = fullfile(figDir, [resDir '_Rest_PDE_FlowMap']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','Flow maps','NumberTitle','off'); imagescn(flow_rest(:,:,:,end), [0 6], [1 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_PDE_Visf']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','PDE Visf','NumberTitle','off'); imagescn(Visf_rest(:,:,:,end), [0 80], [1 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_PDE_PS']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','PDE PS','NumberTitle','off'); imagescn(PS_rest(:,:,:,end), [0 6], [1 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_PDE_E']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','PDE E','NumberTitle','off');; imagescn(E_rest(:,:,:,end), [0 2], [1 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_PDE_Vp']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','PDE Vp','NumberTitle','off');; imagescn(Vp_rest(:,:,:,end), [0 50], [1 slc], scalingFactor); PerfColorMap;
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_Ori']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','Rest Original','NumberTitle','off'); imagescn(ori_rest, [], [1 slc], scalingFactor, 3);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_MOCO']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','Rest MOCO','NumberTitle','off'); imagescn(moco_rest, [], [1 slc], scalingFactor, 3);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_MOCO_NORM']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','Rest MOCO NORM','NumberTitle','off'); imagescn(moco_norm_rest, [], [1 slc], scalingFactor, 3);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_AIF_MOCO']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','AIF MOCO','NumberTitle','off'); imagescn(aif_moco_rest, [], [], 10, 3);
        saveas(h, figName, 'fig')
    end
    
    figName = fullfile(figDir, [resDir '_Rest_PDE_Delay' '.fig']);
    if(onlyReview)
        openfig(figName);
    else
        h = figure('Name','PDE Delay','NumberTitle','off');; imagescn(Delay_rest(:,:,:,end), [0 8], [1 slc], scalingFactor);
        saveas(h, figName, 'fig');
    end
        
    delta = 0.5;

    figName = fullfile(figDir, [resDir '_AIF_Rest_Curves']);

    if(onlyReview)
        openfig(figName);
    else        
        h = figure
        hold on
        plot(delta*[0:numel(aif_rest)-1], aif_rest, 'LineWidth',2)
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
        openfig(figName);
    else        
        h = figure
        hold on
        plot(delta*[0:numel(aif_rest_cin_all_echo0_signal)-1], aif_rest_cin_all_echo0_signal, 'b', 'LineWidth',2);
        plot(delta*[0:numel(aif_rest_cin_all_echo1_signal)-1], aif_rest_cin_all_echo1_signal, 'k', 'LineWidth',2)
        plot(delta*[0:numel(aif_rest_cin_all_echo0_signal_after_R2StarCorrection)-1], aif_rest_cin_all_echo0_signal_after_R2StarCorrection, 'r', 'LineWidth',2)
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
        openfig(figName);
    else        
        
        NN = size(rest_perf, 4);
        
        aif_moco_rest_mask = aif_moco_rest;
        v = repmat(aif_rest_mask, [1 1 NN]);
        ind = find(v>0);
        aif_moco_rest_mask(ind) = 1024;
        
        h = figure('Name','AIF MOCO with Mask, Rest','NumberTitle','off'); imagescn(aif_moco_rest_mask(:,:,1:NN,:), [], [], 10, 3);
        saveas(h, figName, 'fig');
    end
end
