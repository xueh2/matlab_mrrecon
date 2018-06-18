
function [h_flow_stress, h_flow_rest, has_stress, has_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest_WithInfo(resDir, restDir, stressDir, scannerID, patientID, studyID, study_dates_stress, study_time_stress, study_dates_rest, study_time_rest, flow_windowing, onlyReview, checkprocessed, baseDir)
% [h_flow_stress, h_flow_rest, has_stress, has_rest] = PerformGadgetronRecon_Plot_PerfusionCase_StressRest_WithInfo(resDir, restDir, stressDir, scannerID, patientID, studyID, study_dates, flow_windowing, onlyReview, baseDir)
% PerformGadgetronRecon_Plot_PerfusionCase_StressRest_WithInfo('D:\data\perfusion\FAILED', 'D:\data\perfusion\FAILED\meas_MID00761_FID23601_REST_SSFP_Perf_TPAT3_PF34_192x111_new\all_res', 'D:\data\perfusion\FAILED\meas_MID00758_FID23598_STRESS_SSFP_Perf_TPAT3_PF34_192x111_new\all_res', '4096', '1', '100', '20170317')

if(nargin < 8)
    flow_windowing = [0 6];
end

if(nargin < 9)
    onlyReview = 0;
end

if(nargin < 10)
    checkprocessed = 1;
end

if(nargin < 11)
    baseDir = resDir;
end

h_flow_stress = 0;
h_flow_rest = 0;
has_stress = 0;
has_rest = 0;

resDir = [scannerID '_' patientID '_' studyID '_' study_dates_stress];    

ind = find(stressDir=='\');
if(isempty(ind))
    ind = find(stressDir=='/');
end

ind2 = strfind(stressDir, scannerID);
data_name = stressDir(ind(end)+1:ind2(1)-2);
figDir = fullfile(baseDir, study_dates_stress, [data_name '_' resDir '_' study_time_stress '_' study_time_rest '_Figure']);
disp(figDir)
    
if(checkprocessed & ~onlyReview & isFileExist(fullfile(figDir, 'rest.mat')) & isFileExist(fullfile(figDir, 'stress.mat')) & isFileExist(fullfile(figDir, '*Flow*.fig')) )
    disp(['--> Already being processed :' stressDir ' - ' restDir]);
    h_flow_stress = -1;
    h_flow_rest = -1;
    has_stress = 1;
    has_rest = 1;

    return;
end

scalingFactor = 10;

    if(~onlyReview)   
        try
           
            tic
            [rest_perf, ori_rest, moco_rest, moco_norm_rest, PD_rest, input_for_filter_rest, filtered_rest, aif_acq_time_rest, perf_acq_time_rest, dst_acq_time_rest, ... 
                aif_im_rest, aif_moco_rest, aif_rest_cin, aif_rest_cin_Gd, aif_rest_cin_Gd_without_R2Star, aif_rest_cin_Gd_baseline_corrected, ... 
                aif_rest_cin_all_echo0_signal, aif_rest_cin_all_echo1_signal, aif_rest_cin_all_echo0_signal_after_R2StarCorrection, ...
                aif_rest_cin_all_echo0_OverPD_after_R2StarCorrection, aif_rest_cin_all_R2Star,  aif_rest_cin_all_R2Star_SLEP, ... 
                aif_rest_PD, aif_rest_mask, aif_rest_mask_final, aif_rest_LV_mask_plot, aif_rest, aif_rest_baseline_corrected, aif_plots_rest, ... 
                flow_rest, Ki_rest, PS_rest, Vp_rest, Visf_rest, E_rest, SDMap_rest, Delay_rest, ...
                BTEX_Flow_all_rest, BTEX_PS_all_rest, BTEX_Vp_all_rest, BTEX_Visf_all_rest, BTEX_cost_all_rest, ... 
                BTEX_flow_SD_all_rest, BTEX_PS_SD_all_rest, BTEX_Visf_SD_all_rest, BTEX_Vp_SD_all_rest, BTEX_cov_all_rest, ...
                flow_SD_rest, PS_SD_rest, Vp_SD_rest, Visf_SD_rest, BTEX_cov_rest, ...
                CC_F_PS_rest, CC_F_Vp_rest, CC_F_Visf_rest, CC_PS_Vp_rest, CC_PS_Visf_rest, CC_Vp_Visf_rest, ... 
                BTEX_Tc_all_rest, Fermi_Delay_rest] = read_in_GT_Perf_DebugOutput_results(restDir);

            disp(['Load rest - ' num2str(toc)]);
            has_rest = 1;

        catch
            has_rest = 0;
            rest_perf = 0;
            aif_rest = 0;
        end

        try
           
            tic
            [stress_perf, ori_stress, moco_stress, moco_norm_stress, PD_stress, input_for_filter_stress, filtered_stress, aif_acq_time_stress, perf_acq_time_stress, dst_acq_time_stress, ... 
                aif_im_stress, aif_moco_stress, aif_stress_cin, aif_stress_cin_Gd, aif_stress_cin_Gd_without_R2Star, aif_stress_cin_Gd_baseline_corrected, ... 
                aif_stress_cin_all_echo0_signal, aif_stress_cin_all_echo1_signal, aif_stress_cin_all_echo0_signal_after_R2StarCorrection, ...
                aif_stress_cin_all_echo0_OverPD_after_R2StarCorrection, aif_stress_cin_all_R2Star,  aif_stress_cin_all_R2Star_SLEP, ... 
                aif_stress_PD, aif_stress_mask, aif_stress_mask_final, aif_stress_LV_mask_plot, aif_stress, aif_stress_baseline_corrected, aif_plots_stress, ... 
                flow_stress, Ki_stress, PS_stress, Vp_stress, Visf_stress, E_stress, SDMap_stress, Delay_stress, ...
                BTEX_Flow_all_stress, BTEX_PS_all_stress, BTEX_Vp_all_stress, BTEX_Visf_all_stress, BTEX_cost_all_stress, ... 
                BTEX_flow_SD_all_stress, BTEX_PS_SD_all_stress, BTEX_Visf_SD_all_stress, BTEX_Vp_SD_all_stress, BTEX_cov_all_stress, ...
                flow_SD_stress, PS_SD_stress, Vp_SD_stress, Visf_SD_stress, BTEX_cov_stress, ...
                CC_F_PS_stress, CC_F_Vp_stress, CC_F_Visf_stress, CC_PS_Vp_stress, CC_PS_Visf_stress, CC_Vp_Visf_stress, ... 
                BTEX_Tc_all_stress, Fermi_Delay_stress] = read_in_GT_Perf_DebugOutput_results(stressDir);

            disp(['Load stress - ' num2str(toc)]);
            
            has_stress = 1;
        catch
            has_stress = 0;
            stress_perf = 0;
            fb = 0;
            b = 0;
            aif_stress = 0;
        end
    end

    % scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time

    if(~onlyReview)
        if(has_stress || has_rest)  

            if(exist(figDir)==7)
                try
                    rmdir(figDir, 's');
                catch
                end
            end
           
            try
                mkdir(figDir);
            catch
            end
                
            if(has_rest)
                tic
                save(fullfile(figDir, 'rest.mat'), 'rest_perf', 'ori_rest', 'moco_rest', 'moco_norm_rest', 'PD_rest', 'input_for_filter_rest', 'filtered_rest', 'aif_acq_time_rest', 'perf_acq_time_rest', 'dst_acq_time_rest', ... 
                    'aif_im_rest', 'aif_moco_rest', 'aif_rest', 'aif_rest_baseline_corrected', 'aif_rest_cin', ... 
                    'aif_rest_cin_Gd', 'aif_rest_cin_Gd_without_R2Star', 'aif_rest_cin_Gd_baseline_corrected', ... 
                    'aif_rest_cin_all_echo0_signal', 'aif_rest_cin_all_echo1_signal', 'aif_rest_cin_all_echo0_signal_after_R2StarCorrection', ...
                    'aif_rest_cin_all_echo0_OverPD_after_R2StarCorrection', 'aif_rest_cin_all_R2Star', 'aif_rest_cin_all_R2Star_SLEP', 'aif_plots_rest', ...
                    'aif_rest_PD', 'aif_rest_mask', 'aif_rest_mask_final', 'aif_rest_LV_mask_plot', ...
                    'flow_rest', 'Ki_rest', 'PS_rest', 'Vp_rest', 'Visf_rest', 'E_rest', 'SDMap_rest', 'Delay_rest', ...
                    'BTEX_Flow_all_rest', 'BTEX_PS_all_rest', 'BTEX_Visf_all_rest', 'BTEX_Vp_all_rest', 'BTEX_cost_all_rest', ...
                    'BTEX_flow_SD_all_rest', 'BTEX_PS_SD_all_rest', 'BTEX_Visf_SD_all_rest', 'BTEX_Vp_SD_all_rest', 'BTEX_cov_all_rest', ...
                    'flow_SD_rest', 'PS_SD_rest', 'Vp_SD_rest', 'Visf_SD_rest', 'BTEX_cov_rest', ...
                    'CC_F_PS_rest', 'CC_F_Vp_rest', 'CC_F_Visf_rest', 'CC_PS_Vp_rest', 'CC_PS_Visf_rest', 'CC_Vp_Visf_rest', ... 
                    'BTEX_Tc_all_rest', 'Fermi_Delay_rest');
                
                disp(['Save rest - ' num2str(toc)]);
            end

            if(has_stress)
                tic
                save(fullfile(figDir, 'stress.mat'), 'stress_perf', 'ori_stress', 'moco_stress', 'moco_norm_stress', 'PD_stress', 'input_for_filter_stress', 'filtered_stress', 'aif_acq_time_stress', 'perf_acq_time_stress', 'dst_acq_time_stress', ...
                    'aif_im_stress', 'aif_moco_stress', 'aif_stress', 'aif_stress_baseline_corrected', 'aif_stress_cin', ...
                    'aif_stress_cin_Gd', 'aif_stress_cin_Gd_without_R2Star', 'aif_stress_cin_Gd_baseline_corrected', ...
                    'aif_stress_cin_all_echo0_signal', 'aif_stress_cin_all_echo1_signal', 'aif_stress_cin_all_echo0_signal_after_R2StarCorrection', ...
                    'aif_stress_cin_all_echo0_OverPD_after_R2StarCorrection', 'aif_stress_cin_all_R2Star', 'aif_stress_cin_all_R2Star_SLEP', 'aif_plots_stress', ...
                    'aif_stress_PD', 'aif_stress_mask', 'aif_stress_mask_final', 'aif_stress_LV_mask_plot', ...
                    'flow_stress', 'Ki_stress', 'PS_stress', 'Vp_stress', 'Visf_stress', 'E_stress', 'SDMap_stress', 'Delay_stress', ...
                    'BTEX_Flow_all_stress', 'BTEX_PS_all_stress', 'BTEX_Visf_all_stress', 'BTEX_Vp_all_stress', 'BTEX_cost_all_stress', ... 
                    'BTEX_flow_SD_all_stress', 'BTEX_PS_SD_all_stress', 'BTEX_Visf_SD_all_stress', 'BTEX_Vp_SD_all_stress', 'BTEX_cov_all_stress', ...
                    'flow_SD_stress', 'PS_SD_stress', 'Vp_SD_stress', 'Visf_SD_stress', 'BTEX_cov_stress', ...
                    'CC_F_PS_stress', 'CC_F_Vp_stress', 'CC_F_Visf_stress', 'CC_PS_Vp_stress', 'CC_PS_Visf_stress', 'CC_Vp_Visf_stress', ... 
                    'BTEX_Tc_all_stress', 'Fermi_Delay_stress');
                
                disp(['Save stress - ' num2str(toc)]);
            end
            pause(1.0);
        end
    else
        
        tic
        load(fullfile(figDir, 'stress.mat'));
        disp(['Load saved stress - ' num2str(toc)]);
        
        tic
        load(fullfile(figDir, 'rest.mat'));
        disp(['Load saved rest - ' num2str(toc)]);
        
        if(exist(figDir)~=7)
            disp([figDir ' - does not exist']);
            try
                mkdir(figDir);
            catch
            end
        end

        [nameFigs, numFigs] = findFILE(figDir, '*.fig');

        if(numFigs > 0 )
            has_stress = 1;
            has_rest = 1;
        else
            disp([figDir ' - no figures are found']);
            return;
        end
    end

    if(has_stress & has_rest)
    %     a = permute(a, [2 1 3]);
    %     b = permute(b, [2 1 3]);
    % 
    %     h = figure; 
    %     imagescn(cat(4, a, b), [], [], 25)
    %     saveas(h, fullfile(figDir, [resDir '_AIF_FIG']), 'fig')

        if(~onlyReview)
            slc = size(Ki_rest, 4);
            m = size(Ki_rest, 3);

            different_image_size = 0;
            if(size(Ki_rest, 1)~=size(Ki_stress,1) || size(Ki_rest, 2)~=size(Ki_stress,2) )
                disp(['Image size mismatch - Rest :' num2str(size(Ki_rest)) ' - Stress : ' num2str(size(Ki_stress))]);
                different_image_size = 1;
            end

            if(~different_image_size)
                figure; imagescn(cat(4, Ki_stress, Ki_rest), flow_windowing, [m 2*slc], scalingFactor); PerfColorMap;
            end
        else
            slc = 3;
            different_image_size = 0;
        end
        scrsz = get(0, 'ScreenSize');

        C = {'m','g','r', 'k', 'c', 'y', 'b', 'r'};
        
        figName = fullfile(figDir, [resDir '_Stress_Rest_AIF_AcqusitionTimePlot' '.fig']);
        if(onlyReview & isFileExist(figName))

            if(isFileExist(figName))
                openfig(figName);
            end
        else
            h = figure('Name',[ resDir '_AIF Acqusition Time Plot'],'NumberTitle','off', 'Position', [100 100 3*1024 1.5*768]);            
            subplot(2, 1, 1);
            hold on
            
            aif_v = interp1(dst_acq_time_stress, aif_stress_cin_Gd_baseline_corrected, aif_acq_time_stress, 'linear');
            aif_t = (aif_acq_time_stress-dst_acq_time_stress(1))/1000;
            plot(aif_t, aif_v, 'k+', 'MarkerSize',12, 'LineWidth', 2);
            
            perf_acq_time_stress = squeeze(perf_acq_time_stress);
                        
            for ss=1:size(perf_acq_time_stress, 2)
                perf_t = interp1(dst_acq_time_stress, aif_stress_cin_Gd_baseline_corrected, perf_acq_time_stress(:,ss), 'linear');            
                plot((perf_acq_time_stress(:,ss)-dst_acq_time_stress(1))/1000, perf_t, [C{ss} '.'], 'MarkerSize',16);
            end
            xlabel('Acqusition time, in s');
            ylabel('AIF Gd, mmol/L')
            legend('AIF', 'Perf slice 1', 'Perf slice 2', 'Perf slice 3')
            title('Stress')
            
            xlim_v = [min(aif_t)-0.5 max(aif_t)+0.5];
            ylim_v = [min(aif_stress_cin_Gd_baseline_corrected)-1 max(aif_stress_cin_Gd_baseline_corrected)+1];
            
            plot((dst_acq_time_stress-dst_acq_time_stress(1))/1000, aif_stress_cin_Gd_baseline_corrected, 'b--');            

            plot([xlim_v(1) xlim_v(2)], [ylim_v(1)+0.5 ylim_v(1)+0.5], 'b-')
            plot(aif_t, ylim_v(1)+0.5, 'k+', 'MarkerSize',12, 'LineWidth', 2)
            
            xlim(xlim_v);
            ylim(ylim_v);
            hold off
            box on
            grid on
            
            % ------------------------
            
            subplot(2, 1, 2);
            hold on
            aif_v = interp1(dst_acq_time_rest, aif_rest_cin_Gd_baseline_corrected, aif_acq_time_rest, 'linear');
            aif_t = (aif_acq_time_rest-dst_acq_time_rest(1))/1000;
            plot(aif_t, aif_v, 'k+', 'MarkerSize',12, 'LineWidth', 2);
            
            perf_acq_time_rest = squeeze(perf_acq_time_rest);
            
            for ss=1:size(perf_acq_time_rest, 2)
                perf_t = interp1(dst_acq_time_rest, aif_rest_cin_Gd_baseline_corrected, perf_acq_time_rest(:,ss), 'linear');            
                plot((perf_acq_time_rest(:,ss)-dst_acq_time_rest(1))/1000, perf_t, [C{ss} '.'], 'MarkerSize',16);
            end
            xlabel('Acqusition time, in s');
            ylabel('AIF Gd, mmol/L')
            legend('AIF', 'Perf slice 1', 'Perf slice 2', 'Perf slice 3')
            title('Rest')
            
            xlim_v = [min(aif_t)-0.5 max(aif_t)+0.5];
            ylim_v = [min(aif_rest_cin_Gd_baseline_corrected)-1 max(aif_rest_cin_Gd_baseline_corrected)+1];
            
            plot((dst_acq_time_rest-dst_acq_time_rest(1))/1000, aif_rest_cin_Gd_baseline_corrected, 'b--');            

            plot([xlim_v(1) xlim_v(2)], [ylim_v(1)+0.5 ylim_v(1)+0.5], 'b-')
            plot(aif_t, ylim_v(1)+0.5, 'k+', 'MarkerSize',12, 'LineWidth', 2)
            
            xlim(xlim_v);
            ylim(ylim_v);            
            hold off
            box on
            grid on
            
            saveas(h, figName, 'fig');
        end            
        

        figName = fullfile(figDir, [resDir '_Stress_Rest_Ki' '.fig']);
        if(onlyReview & isFileExist(figName))

            if(isFileExist(figName))
                openfig(figName);

                for s=1:slc
                    figName = fullfile(figDir, [resDir '_PDE_StressKiMap - ' num2str(s)]);
                    openfig(figName);
                end

                for s=1:slc
                    figName = fullfile(figDir, [resDir '_PDE_RestKiMap - ' num2str(s)]);
                    openfig(figName);
                end
            end
        else
            numKi = size(Ki_stress, 3);
            
            if(~different_image_size)
                h = figure('Name',[ resDir '_Ki Maps'],'NumberTitle','off'); imagescn(cat(4, Ki_stress, Ki_rest), flow_windowing, [numKi 2*slc], scalingFactor); PerfColorMap;
                saveas(h, figName, 'fig');
            end

            for s=1:slc
                h = figure('Name', [ resDir '_Stress Ki maps - ' num2str(s)],'NumberTitle','off'); imagescn(Ki_stress(:,:,:,s), flow_windowing, [1 numKi], scalingFactor); PerfColorMap;
                figName = fullfile(figDir, [resDir '_PDE_StressKiMap - ' num2str(s) '.fig']);
                saveas(h, figName, 'fig');
            end

            for s=1:slc
                h = figure('Name',[ resDir '_Rest Ki maps - ' num2str(s)],'NumberTitle','off'); imagescn(Ki_rest(:,:,:,s), flow_windowing, [1 numKi], scalingFactor); PerfColorMap;
                figName = fullfile(figDir, [resDir '_PDE_RestKiMap - ' num2str(s) '.fig']);
                saveas(h, figName, 'fig');
            end
        end

        figName = fullfile(figDir, [resDir '_Stress_Rest_Gd' '.fig']);
        if(onlyReview & isFileExist(figName))
            % openfig(figName);
        else
            NS = size(stress_perf, 4);
            NR = size(rest_perf, 4);

            NN = NS;
            if(NN>NR) NN = NR; end

            if(~different_image_size)
                h = figure('Name',[ resDir '_Rest-Stress Gd'],'NumberTitle','off'); imagescn(cat(3, stress_perf(:,:,:,1:NN), rest_perf(:,:,:,1:NN)), [0 1.5], [2 slc], 10, 4);
                saveas(h, figName, 'fig');
            else
                h = figure('Name',[ resDir '_Rest Gd'],'NumberTitle','off'); imagescn(rest_perf(:,:,:,1:NN), [0 1.5], [1 slc], 10, 4);
                saveas(h, fullfile(figDir, [resDir '_Rest_Gd' '.fig']), 'fig');

                h = figure('Name',[ resDir '_Stress Gd'],'NumberTitle','off'); imagescn(stress_perf(:,:,:,1:NN), [0 1.5], [1 slc], 10, 4);
                saveas(h, fullfile(figDir, [resDir '_Stress_Gd' '.fig']), 'fig');
            end
        end

        figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_FlowMap' '.fig']);
        if(onlyReview & isFileExist(figName))
            if(isFileExist(figName)) 
                openfig(figName); 
            end
        else
            if(~different_image_size)
                h = figure('Name',[ resDir '_Flow maps'],'NumberTitle','off'); imagescn(cat(3, flow_stress, flow_rest), flow_windowing, [2 slc], scalingFactor); PerfColorMap;
                saveas(h, figName, 'fig');
            end
        end

        figName = fullfile(figDir, [resDir '_PDE_StressFlowMap' '.fig']);
        figName1 = fullfile(figDir, [resDir '_PDE_StressFlowMap - ' num2str(1)]);
        figName2 = fullfile(figDir, [resDir '_PDE_StressFlowMap - ' num2str(2)]);
        figName3 = fullfile(figDir, [resDir '_PDE_StressFlowMap - ' num2str(3)]);

        if(onlyReview & isFileExist(figName))
            openfig(figName);
            h_flow_stress(1) = openfig(figName1);
            h_flow_stress(2) = openfig(figName2);
            h_flow_stress(3) = openfig(figName3);
        else
            h = figure('Name',[ resDir '_Stress Flow maps'],'NumberTitle','off'); imagescn(flow_stress(:,:,:,end), flow_windowing, [1 slc], scalingFactor); PerfColorMap;
            saveas(h, figName, 'fig');

            for s=1:slc
                h_flow_stress(s) = figure('Name', [ resDir '_Stress Flow maps - ' num2str(s)],'NumberTitle','off'); imagescn(flow_stress(:,:,s,end), flow_windowing, [1 1], scalingFactor); PerfColorMap;
                figName = fullfile(figDir, [resDir '_PDE_StressFlowMap - ' num2str(s) '.fig']);
                saveas(h_flow_stress(s), figName, 'fig');
            end
        end

        figName = fullfile(figDir, [resDir '_PDE_RestFlowMap' '.fig']);
        figName1 = fullfile(figDir, [resDir '_PDE_RestFlowMap - ' num2str(1) '.fig']);
        figName2 = fullfile(figDir, [resDir '_PDE_RestFlowMap - ' num2str(2) '.fig']);
        figName3 = fullfile(figDir, [resDir '_PDE_RestFlowMap - ' num2str(3) '.fig']);
        if(onlyReview & isFileExist(figName))
            openfig(figName);
            h_flow_rest(1) = openfig(figName1);
            h_flow_rest(2) = openfig(figName2);
            h_flow_rest(3) = openfig(figName3);
        else
            h = figure('Name',[ resDir '_Rest Flow maps'],'NumberTitle','off'); imagescn(flow_rest(:,:,:,end), flow_windowing, [1 slc], scalingFactor); PerfColorMap;
            saveas(h, figName, 'fig');

            for s=1:slc
                h_flow_rest(s) = figure('Name',[ resDir '_Rest Flow maps - ' num2str(s)],'NumberTitle','off'); imagescn(flow_rest(:,:,s,end), flow_windowing, [1 1], scalingFactor); PerfColorMap;
                figName = fullfile(figDir, [resDir '_PDE_RestFlowMap - ' num2str(s) '.fig']);
                saveas(h_flow_rest(s), figName, 'fig');
            end
        end

        figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_Visf' '.fig']);   
        % if(~isFileExist(figName) & ~different_image_size)
        if(~different_image_size)

            if(onlyReview & isFileExist(figName))
                openfig(figName);
            else        
                h = figure('Name',[ resDir '_PDE Visf'],'NumberTitle','off'); imagescn(cat(3, Visf_stress(:,:,:,end), Visf_rest(:,:,:,end)), [0 80], [2 slc], scalingFactor); ECVColorMap;
                saveas(h, figName, 'fig');
            end

            figName = fullfile(figDir, [resDir '_PDE_Stress_Visf' '.fig']);
            figName1 = fullfile(figDir, [resDir '_PDE_Stress_Visf_Map - ' num2str(1)  '.fig']);
            figName2 = fullfile(figDir, [resDir '_PDE_Stress_Visf_Map - ' num2str(2)  '.fig']);
            figName3 = fullfile(figDir, [resDir '_PDE_Stress_Visf_Map - ' num2str(3)  '.fig']);

            if(onlyReview & isFileExist(figName) & isFileExist(figName1))
                openfig(figName);
                h_visf_stress(1) = openfig(figName1);
                h_visf_stress(2) = openfig(figName2);
                h_visf_stress(3) = openfig(figName3);
            else
                h = figure('Name',[ resDir '_Stress Visf maps'],'NumberTitle','off'); imagescn(Visf_stress(:,:,:,end), [0 80], [1 slc], scalingFactor); ECVColorMap;
                saveas(h, figName, 'fig');

                for s=1:slc
                    h_visf_stress(s) = figure('Name', [ resDir '_Stress Visf maps - ' num2str(s)],'NumberTitle','off'); imagescn(Visf_stress(:,:,s,end), [0 80], [1 1], scalingFactor); ECVColorMap;
                    figName = fullfile(figDir, [resDir '_PDE_Stress_Visf_Map - ' num2str(s) '.fig']);
                    saveas(h_visf_stress(s), figName, 'fig');
                end
            end
            
            figName = fullfile(figDir, [resDir '_PDE_Rest_Visf' '.fig']);
            figName1 = fullfile(figDir, [resDir '_PDE_Rest_Visf_Map - ' num2str(1)  '.fig']);
            figName2 = fullfile(figDir, [resDir '_PDE_Rest_Visf_Map - ' num2str(2)  '.fig']);
            figName3 = fullfile(figDir, [resDir '_PDE_Rest_Visf_Map - ' num2str(3)  '.fig']);

            if(onlyReview & isFileExist(figName) & isFileExist(figName1) )
                openfig(figName);
                h_visf_rest(1) = openfig(figName1);
                h_visf_rest(2) = openfig(figName2);
                h_visf_rest(3) = openfig(figName3);
            else
                h = figure('Name',[ resDir '_Rest Visf maps'],'NumberTitle','off'); imagescn(Visf_rest(:,:,:,end), [0 80], [1 slc], scalingFactor); ECVColorMap;
                saveas(h, figName, 'fig');

                for s=1:slc
                    h_visf_rest(s) = figure('Name', [ resDir '_Rest Visf maps - ' num2str(s)],'NumberTitle','off'); imagescn(Visf_rest(:,:,s,end), [0 80], [1 1], scalingFactor); ECVColorMap;
                    figName = fullfile(figDir, [resDir '_PDE_Rest_Visf_Map - ' num2str(s) '.fig']);
                    saveas(h_visf_rest(s), figName, 'fig');
                end
            end
            
            
            
            figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_PS' '.fig']);
            if(onlyReview & isFileExist(figName))
                openfig(figName);
            else
                h = figure('Name',[ resDir '_PDE PS'],'NumberTitle','off'); imagescn(cat(3, PS_stress(:,:,:,end), PS_rest(:,:,:,end)), [0 2], [2 slc], scalingFactor); PSColorMap;
                saveas(h, figName, 'fig');
            end

            figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_E' '.fig']);
            if(onlyReview & isFileExist(figName))
                openfig(figName);
            else
                h = figure('Name',[ resDir '_PDE E'],'NumberTitle','off');; imagescn(cat(3, E_stress(:,:,:,end), E_rest(:,:,:,end)), [0 2], [2 slc], scalingFactor); PerfColorMap;
                saveas(h, figName, 'fig');
            end

            figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_Vb' '.fig']);
            if(onlyReview & isFileExist(figName))
                openfig(figName);
            else
                h = figure('Name',[ resDir '_PDE Vb'],'NumberTitle','off');; imagescn(cat(3, Vp_stress(:,:,:,end), Vp_rest(:,:,:,end)), [0 20], [2 slc], scalingFactor); MBVColorMap;
                saveas(h, figName, 'fig');
            end

            figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_Tc' '.fig']);
            if(onlyReview & isFileExist(figName))
                openfig(figName);
            else
                h = figure('Name',[ resDir '_PDE Tc'],'NumberTitle','off');; imagescn(cat(3, BTEX_Tc_all_stress, BTEX_Tc_all_rest), [0 10], [2 slc], scalingFactor); PerfColorMap;
                saveas(h, figName, 'fig');
            end
            
            try
                figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_Flow_SD' '.fig']);
                if(onlyReview & isFileExist(figName))
                    openfig(figName);
                else
                    h = figure('Name',[ resDir '_PDE Flow SD'],'NumberTitle','off');; imagescn(cat(3, flow_SD_stress(:,:,:,end), flow_SD_rest(:,:,:,end)), [0 1], [2 slc], scalingFactor); PerfColorMap;
                    saveas(h, figName, 'fig');
                end
            catch
            end
            
            try
                figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_PS_SD' '.fig']);
                if(onlyReview & isFileExist(figName))
                    openfig(figName);
                else
                    h = figure('Name',[ resDir '_PDE PS SD'],'NumberTitle','off');; imagescn(cat(3, PS_SD_stress(:,:,:,end), PS_SD_rest(:,:,:,end)), [0 1], [2 slc], scalingFactor); PSColorMap;
                    saveas(h, figName, 'fig');
                end
            catch
            end
            
            try
                figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_Visf_SD' '.fig']);
                if(onlyReview & isFileExist(figName))
                    openfig(figName);
                else
                    h = figure('Name',[ resDir '_PDE Visf SD'],'NumberTitle','off');; imagescn(cat(3, 100*Visf_SD_stress(:,:,:,end), 100*Visf_SD_rest(:,:,:,end)), [0 10], [2 slc], scalingFactor); ECVColorMap;
                    saveas(h, figName, 'fig');
                end
            catch
            end
            
            try
                figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_Vp_SD' '.fig']);
                if(onlyReview & isFileExist(figName))
                    openfig(figName);
                else
                    h = figure('Name',[ resDir '_PDE Vp SD'],'NumberTitle','off');; imagescn(cat(3, 100*Vp_SD_stress(:,:,:,end), 100*Vp_SD_rest(:,:,:,end)), [0 4], [2 slc], scalingFactor); MBVColorMap;
                    saveas(h, figName, 'fig');
                end
            catch
            end
            
            try                  
                for slc=1:SLC
                    h = figure('Name',[ resDir '_PDE Flow to PS/Vp/Visf correlation coefficients' '_SLC' num2str(slc)],'NumberTitle','off'); imagescn(cat(3, CC_F_PS_stress(:,:,slc), CC_F_Vp_stress(:,:,slc), CC_F_Visf_stress(:,:,slc)), [-1 1], [1 3], scalingFactor); PerfColorMap;
                    saveas(h, figName, 'fig');

                    h = figure('Name',[ resDir '_PDE PS to Flow/Vp/Visf correlation coefficients' '_SLC' num2str(slc)],'NumberTitle','off'); imagescn(cat(3, CC_F_PS_stress(:,:,slc), CC_PS_Vp_stress(:,:,slc), CC_PS_Visf_stress(:,:,slc)), [-1 1], [1 3], scalingFactor); PerfColorMap;
                    saveas(h, figName, 'fig');

                    h = figure('Name',[ resDir '_PDE Vp to Flow/PS/Visf correlation coefficients' '_SLC' num2str(slc)],'NumberTitle','off'); imagescn(cat(3, CC_F_Vp_stress(:,:,slc), CC_PS_Vp_stress(:,:,slc), CC_Vp_Visf_stress(:,:,slc)), [-1 1], [1 3], scalingFactor); PerfColorMap;
                    saveas(h, figName, 'fig');

                    h = figure('Name',[ resDir '_PDE Visf to Flow/PS/Vp correlation coefficients' '_SLC' num2str(slc)],'NumberTitle','off'); imagescn(cat(3, CC_F_Visf_stress(:,:,slc), CC_PS_Visf_stress(:,:,slc), CC_Vp_Visf_stress(:,:,slc)), [-1 1], [1 3], scalingFactor); PerfColorMap;
                    saveas(h, figName, 'fig');
                end
            catch
            end
            
            figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_Delay' '.fig']);
            if(onlyReview & isFileExist(figName))
                openfig(figName);
            else
                h = figure('Name',[ resDir '_PDE Delay'],'NumberTitle','off');; imagescn(cat(3, Delay_stress(:,:,:,end), Delay_rest(:,:,:,end)), [0 8], [2 slc], scalingFactor); PerfColorMap;
                saveas(h, figName, 'fig');
            end

            if(~isempty(Fermi_Delay_stress) & ~isempty(Fermi_Delay_rest))
                figName = fullfile(figDir, [resDir '_Stress_Rest_PDE_Fermi_Delay' '.fig']);
                if(onlyReview & isFileExist(figName))
                    openfig(figName);
                else
                    h = figure('Name',[ resDir '_Fermi Delay'],'NumberTitle','off');; imagescn(cat(3, Fermi_Delay_stress(:,:,:,end), Fermi_Delay_rest(:,:,:,end)), [0 8], [2 slc], scalingFactor); PerfColorMap;
                    saveas(h, figName, 'fig');
                end
            end
            
            if(onlyReview)
                return;
            end
            
            NS = size(ori_stress, 3);
            NR = size(ori_stress, 3);

            NN = NS;
            if(NN>NR) NN = NR; end
            
            figName = fullfile(figDir, [resDir '_Stress_Rest_Ori' '.fig']);
            if(onlyReview & isFileExist(figName))
                openfig(figName);
            else
                h = figure('Name',[ resDir '_Rest-Stress Original'],'NumberTitle','off'); imagescn(cat(4, ori_stress(:,:,1:NN,:), ori_rest(:,:,1:NN,:)), [], [2 slc], scalingFactor, 3);
                saveas(h, figName, 'fig');
            end

            figName = fullfile(figDir, [resDir '_Stress_Rest_MOCO' '.fig']);
            if(onlyReview & isFileExist(figName))
                if(isFileExist(figName) ) 
                    openfig(figName);
                else
                    h = figure('Name',[ resDir '_Rest-Stress MOCO'],'NumberTitle','off'); imagescn(cat(4, moco_stress(:,:,1:NN,:), moco_rest(:,:,1:NN,:)), [], [2 slc], scalingFactor, 3);
                    saveas(h, figName, 'fig');
                end
            else
                h = figure('Name',[ resDir '_Rest-Stress MOCO'],'NumberTitle','off'); imagescn(cat(4, moco_stress(:,:,1:NN,:), moco_rest(:,:,1:NN,:)), [], [2 slc], scalingFactor, 3);
                saveas(h, figName, 'fig');
            end

            figName = fullfile(figDir, [resDir '_Stress_Rest_CoMOCO' '.fig']);
            if(onlyReview & isFileExist(figName))
                if(isFileExist(figName) ) 
                    openfig(figName);
                else
                    ss = moco_stress(:,:,[1 4 NN],:);
                    rr = moco_rest(:,:,[1 4 NN],:);

                    h = figure('Name',[ resDir '_Rest-Stress SR-PD CoMOCO'],'NumberTitle','off'); imagescn(cat(4, ss, rr), [], [2*slc 3], scalingFactor);
                    saveas(h, figName, 'fig');
                end
            else
                ss = moco_stress(:,:,[1 4 NN],:);
                rr = moco_rest(:,:,[1 4 NN],:);

                h = figure('Name',[ resDir '_Rest-Stress SR-PD CoMOCO'],'NumberTitle','off'); imagescn(cat(4, ss, rr), [], [2*slc 3], scalingFactor);
                saveas(h, figName, 'fig');
            end

            figName = fullfile(figDir, [resDir '_Stress_Rest_MOCO_NORM' '.fig']);
            if(onlyReview & isFileExist(figName))
                if(isFileExist(figName) ) 
                    openfig(figName);
                else
                    h = figure('Name',[ resDir '_Rest-Stress MOCO NORM'],'NumberTitle','off'); imagescn(cat(4, moco_norm_stress(:,:,1:NN,:), moco_norm_rest(:,:,1:NN,:)), [0 3], [2 slc], scalingFactor, 3);
                    saveas(h, figName, 'fig');
                end
            else
                h = figure('Name',[ resDir '_Rest-Stress MOCO NORM'],'NumberTitle','off'); imagescn(cat(4, moco_norm_stress(:,:,1:NN,:), moco_norm_rest(:,:,1:NN,:)), [0 3], [2 slc], scalingFactor, 3);
                saveas(h, figName, 'fig');
            end
        end

        if(~isempty(input_for_filter_stress))
            figName = fullfile(figDir, [resDir '_Stress_MOCO_FIL' '.fig']);
            if(onlyReview & isFileExist(figName))
                if(isFileExist(figName) ) 
                    openfig(figName);
                else
                    h = figure('Name',[ resDir '_Stress MOCO with/without filtering'],'NumberTitle','off'); imagescn(cat(4, input_for_filter_stress(:,:,1:NN,:), filtered_stress(:,:,1:NN,:)), [], [2 slc], scalingFactor, 3);
                    saveas(h, figName, 'fig');
                end
            else
                h = figure('Name',[ resDir '_Stress MOCO with/without filtering'],'NumberTitle','off'); imagescn(cat(4, input_for_filter_stress(:,:,1:NN,:), filtered_stress(:,:,1:NN,:)), [], [2 slc], scalingFactor, 3);
                saveas(h, figName, 'fig');
            end
        end
        
        if(~isempty(input_for_filter_rest))
            figName = fullfile(figDir, [resDir '_Rest_MOCO_FIL' '.fig']);
            if(onlyReview & isFileExist(figName))
                if(isFileExist(figName) ) 
                    openfig(figName);
                else
                    h = figure('Name',[ resDir '_Rest MOCO with/without filtering'],'NumberTitle','off'); imagescn(cat(4, input_for_filter_rest(:,:,1:NN,:), filtered_rest(:,:,1:NN,:)), [], [2 slc], scalingFactor, 3);
                    saveas(h, figName, 'fig');
                end
            else
                h = figure('Name',[ resDir '_Rest MOCO with/without filtering'],'NumberTitle','off'); imagescn(cat(4, input_for_filter_rest(:,:,1:NN,:), filtered_rest(:,:,1:NN,:)), [], [2 slc], scalingFactor, 3);
                saveas(h, figName, 'fig');
            end
        end
        
        if(~different_image_size)
            figName = fullfile(figDir, [resDir '_Stress_Rest_AIF_ORI' '.fig']);
            if(onlyReview & isFileExist(figName))
                openfig(figName);
            else        
                h = figure('Name',[ resDir '_AIF Original'],'NumberTitle','off'); imagescn(cat(3, aif_im_stress, aif_im_rest), [], [], 10, 4);
                saveas(h, figName, 'fig');
            end
        end

        if(~different_image_size)
            figName = fullfile(figDir, [resDir '_Stress_Rest_AIF_MOCO' '.fig']);
            if(onlyReview & isFileExist(figName))
                openfig(figName);
            else        
                h = figure('Name',[ resDir '_AIF MOCO'],'NumberTitle','off'); imagescn(cat(4, aif_moco_stress(:,:,1:NN,:), aif_moco_rest(:,:,1:NN,:)), [], [], 10, 3);
                saveas(h, figName, 'fig');
            end
        end

        figName = fullfile(figDir, [resDir '_Stress_AIF_Mask' '.fig']);
        if(onlyReview & isFileExist(figName))
            openfig(figName);
        else        

            aif_moco_stress_mask = aif_moco_stress;
            v = repmat(aif_stress_mask, [1 1 NN]);
            ind = find(v>0);
            aif_moco_stress_mask(ind) = 1024;

            h = figure('Name',[ resDir '_AIF MOCO with Mask, Stress'],'NumberTitle','off'); imagescn(aif_moco_stress_mask(:,:,1:NN,:), [], [], 10, 3);
            saveas(h, figName, 'fig');
        end

        figName = fullfile(figDir, [resDir '_Rest_AIF_Mask' '.fig']);
        if(onlyReview & isFileExist(figName))
            openfig(figName);
        else        

            aif_moco_rest_mask = aif_moco_rest;
            v = repmat(aif_rest_mask, [1 1 NN]);
            ind = find(v>0);
            aif_moco_rest_mask(ind) = 1024;

            h = figure('Name',[ resDir '_AIF MOCO with Mask, Rest'],'NumberTitle','off'); imagescn(aif_moco_rest_mask(:,:,1:NN,:), [], [], 10, 3);
            saveas(h, figName, 'fig');
        end

        % AIF mask on perf
        figName = fullfile(figDir, [resDir '_Stress_With_AIF_Mask_Symbol' '.fig']);
        figName2 = fullfile(figDir, [resDir '_Rest_With_AIF_Mask_Symbol' '.fig']);
        if(onlyReview & isFileExist(figName))
            openfig(figName);
            openfig(figName2);
        else
            h = figure('Name',[ resDir '_AIF mask, Stress'],'NumberTitle','off'); imagescn(aif_stress_LV_mask_plot, [], [], 10);
            saveas(h, figName, 'fig');
            
            h = figure('Name',[ resDir '_AIF mask, Rest'],'NumberTitle','off'); imagescn(aif_rest_LV_mask_plot, [], [], 10);
            saveas(h, figName2, 'fig');
        end
        
%         figName = fullfile(figDir, [resDir '_Stress_MOCO_With_AIF_Mask' '.fig']);
%         % if(onlyReview & isFileExist(figName))
%         if(0)
%             openfig(figName);
%         else        
% 
%             [indx, indy] = find(aif_stress_mask>0);
%             c = round([mean(indx) mean(indy)]);
%             
%             maxR = 0;
%             for ii=1:numel(indx)                
%                 r = norm([indx(ii)-c(1) indy(ii)-c(2)]);
%                 if(r>maxR)
%                     maxR = r;
%                 end                
%             end
%             
%             if(maxR<4) maxR = 4; end
%             
%             RO = size(moco_norm_stress, 1);
%             E1 = size(moco_norm_stress, 2);
%             
%             RO_ratio = RO / size(aif_stress_mask, 1);
%             E1_ratio = E1 / size(aif_stress_mask, 2);
%             
%             mask = zeros(RO, E1);
%             
%             c_perf = [c(1)*RO_ratio c(2)*E1_ratio];            
%             maxR_perf = maxR * max([RO_ratio E1_ratio]);
%                         
%             R_times = 4;
%             
%             rx = maxR_perf*R_times;
%             ry = rx*E1/RO;
% 
%             for e1=1:E1
%                 for ro=1:RO   
%                     dx = abs(ro-c_perf(1));
%                     dy = abs(e1-c_perf(2));
%                     d = norm([dx dy]);                   
% 
%                     tt = dx*dx/(rx*rx) + dy*dy/(ry*ry);
%                     
%                     if(tt<=1)
%                         mask(ro, e1) = 1;
%                     end
%                 end
%             end
%            
%             mask_boundary = zeros(RO, E1);
%             for e1=2:E1-1
%                 for ro=2:RO-1
%                     a = mask(ro-1:ro+1, e1-1:e1+1);
%                     ind = find(a==0);
%                     ind2 = find(a>0);
%                     
%                     if(numel(ind)>0 & numel(ind2)>0)
%                         mask_boundary(ro, e1) = 1;
%                     end
%                 end
%             end
%             
%             SLC = size(moco_norm_stress, 4);
%             N = size(moco_norm_stress, 3);
%             
%             mI = max(moco_norm_stress(:));
% 
%             moco_norm = moco_norm_stress;
%             
%             for e1=1:E1
%                 for ro=1:RO                    
%                     if(mask_boundary(ro, e1)==1)
%                         for slc=1:SLC
%                             for n=1:N
%                                 moco_norm(ro, e1, n, slc) = mI+1;
%                             end
%                         end
%                     end
%                 end
%             end
%             
%             h = figure('Name',[ resDir '_Stress MOCO with AIF mask, Stress'],'NumberTitle','off'); imagescn(moco_norm, [], [], 10, 3);
%             saveas(h, figName, 'fig');
%         end
%         
%         figName = fullfile(figDir, [resDir '_Rest_MOCO_With_AIF_Mask' '.fig']);
%         %if(onlyReview & isFileExist(figName))
%         if(0)
%             openfig(figName);
%         else        
% 
%             [indx, indy] = find(aif_rest_mask>0);
%             c = round([mean(indx) mean(indy)]);
%             
%             maxR = 0;
%             for ii=1:numel(indx)                
%                 r = norm([indx(ii)-c(1) indy(ii)-c(2)]);
%                 if(r>maxR)
%                     maxR = r;
%                 end                
%             end
%             
%             RO = size(moco_norm_rest, 1);
%             E1 = size(moco_norm_rest, 2);
%             
%             RO_ratio = RO / size(aif_rest_mask, 1);
%             E1_ratio = E1 / size(aif_rest_mask, 2);
%             
%             mask = zeros(RO, E1);
%             
%             c_perf = [c(1)*RO_ratio c(2)*E1_ratio];            
%             maxR_perf = maxR * max([RO_ratio E1_ratio]);
%             
%             R_times = 4;
%             
%             rx = maxR_perf*R_times;
%             ry = rx*E1/RO;
% 
%             for e1=1:E1
%                 for ro=1:RO   
%                     dx = abs(ro-c_perf(1));
%                     dy = abs(e1-c_perf(2));
%                     d = norm([dx dy]);                   
% 
%                     tt = dx*dx/(rx*rx) + dy*dy/(ry*ry);
%                     
%                     if(tt<=1)
%                         mask(ro, e1) = 1;
%                     end
%                 end
%             end
%            
%             mask_boundary = zeros(RO, E1);
%             for e1=2:E1-1
%                 for ro=2:RO-1
%                     a = mask(ro-1:ro+1, e1-1:e1+1);
%                     ind = find(a==0);
%                     ind2 = find(a>0);
%                     
%                     if(numel(ind)>0 & numel(ind2)>0)
%                         mask_boundary(ro, e1) = 1;
%                     end
%                 end
%             end
%             
%             SLC = size(moco_norm_rest, 4);
%             N = size(moco_norm_rest, 3);
%             
%             mI = max(moco_norm_rest(:));
% 
%             moco_norm = moco_norm_rest;
%             
%             for e1=1:E1
%                 for ro=1:RO                    
%                     if(mask_boundary(ro, e1)==1)
%                         for slc=1:SLC
%                             for n=1:N
%                                 moco_norm(ro, e1, n, slc) = mI+1;
%                             end
%                         end
%                     end
%                 end
%             end
%             
%             h = figure('Name',[ resDir '_Rest MOCO with AIF mask, Rest'],'NumberTitle','off'); imagescn(moco_norm, [], [], 10, 3);
%             saveas(h, figName, 'fig');
%         end
        
        %% other line plots
        delta = 0.5;

        figName = fullfile(figDir, [resDir '_AIF_Stress_Rest_Curves' '.fig']);
        if(onlyReview & isFileExist(figName))
            openfig(figName);
        else        
            h = figure('Name',[ resDir '_AIF_Stress_Rest_Curves'],'NumberTitle','off');
            hold on
            plot(delta*[0:numel(aif_stress)-1], aif_stress, 'r', 'LineWidth',2);
            plot(delta*[0:numel(aif_rest)-1], aif_rest, 'LineWidth',2)
            hold off
            set(gcf, 'Units', 'normalized', 'Position', [0.2 0.2 0.5 0.5])
            box on
            grid on
            grid MINOR
            legend('stress', 'rest')
            xlabel('second')
            ylabel('Gd [mmol/ml]')
            title('AIF in Gd')

            saveas(h, figName, 'fig');
        end

        figName = fullfile(figDir, [resDir '_AIF_Stress_Curves' '.fig']);
        if(onlyReview & isFileExist(figName))
            openfig(figName);
        else        
            h = figure('Name',[ resDir '_AIF_Stress_Curves'],'NumberTitle','off');
            hold on
            plot(delta*[0:numel(aif_stress_cin_all_echo0_signal)-1], aif_stress_cin_all_echo0_signal, 'b', 'LineWidth',2);
            plot(delta*[0:numel(aif_stress_cin_all_echo1_signal)-1], aif_stress_cin_all_echo1_signal, 'k', 'LineWidth',2)
            plot(delta*[0:numel(aif_stress_cin_all_echo0_signal_after_R2StarCorrection)-1], aif_stress_cin_all_echo0_signal_after_R2StarCorrection, 'r', 'LineWidth',2)
            hold off
            set(gcf, 'Units', 'normalized', 'Position', [0.2 0.2 0.5 0.5])
            box on
            grid on
            grid MINOR
            legend('first echo', 'second echo', 'first echo with T2* correction')
            xlabel('second')
            ylabel('internsity')
            title('AIF T2* correction, stress')

            saveas(h, figName, 'fig');
        end

        figName = fullfile(figDir, [resDir '_AIF_Rest_Curves' '.fig']);
        if(onlyReview & isFileExist(figName))
            openfig(figName);
        else        
            h = figure('Name',[ resDir '_AIF_Rest_Curves'],'NumberTitle','off');
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
        
        figName = fullfile(figDir, [resDir '_AIF_Stress_Rest_Plots' '.fig']);
        if(onlyReview & isFileExist(figName))
            openfig(figName);
        else        
            h = figure('Name',[ resDir '_AIF_Stress_Rest_Plots'],'NumberTitle','off');
            imagescn(cat(4, aif_plots_stress, aif_plots_rest), [], [2 3], [26], 1);
            saveas(h, figName, 'fig');
        end
    else
        if(has_stress & ~has_rest)
    %         b = permute(b, [2 1 3]);
    % 
    %         h = figure; 
    %         imagescn(b, [], [], 25)
    %         saveas(h, fullfile(figDir, [resDir '_AIF_Stress_FIG']), 'fig')

            slc = size(flow_stress, 3);
            m = size(flow_stress, 4);

            figure; imagescn(flow_stress, flow_windowing, [m slc], scalingFactor); PerfColorMap;

            figure; imagescn(stress_perf, [0 3], [], scalingFactor, 4);

            delta = 0.5;

            if(~isempty(aif_stress))
                figure
                hold on
                plot(delta*[0:numel(aif_stress)-1], aif_stress, 'r', 'LineWidth',2);
                hold off
                set(gcf, 'Units', 'normalized', 'Position', [0.2 0.2 0.5 0.5])
                box on
                grid on
                grid MINOR
                legend('stress')
                xlabel('second')
                ylabel('Gd [mmol/ml]')
                title('AIF in Gd')
            end
        else
            if(~has_stress & has_rest)
    %             a = permute(a, [2 1 3]);
    % 
    %             h = figure; 
    %             imagescn(a, [], [], 25)
    %             saveas(h, fullfile(figDir, [resDir '_AIF_Rest_FIG']), 'fig')

%                 slc = size(fa, 3);
%                 m = size(fa, 4);
% 
%                 figure; imagescn(fa, flow_windowing, [m slc], scalingFactor); PerfColorMap;

                figure; imagescn(rest_perf, [0 3], [], 10, 4);

                delta = 0.5;

                if(~isempty(aif_rest))
                    figure
                    hold on
                    plot(delta*[0:numel(aif_rest)-1], aif_rest, 'r', 'LineWidth',2);
                    hold off
                    set(gcf, 'Units', 'normalized', 'Position', [0.2 0.2 0.5 0.5])
                    box on
                    grid on
                    grid MINOR
                    legend('rest')
                    xlabel('second')
                    ylabel('Gd [mmol/ml]')
                    title('AIF in Gd')
                end
            end
        end
    end
end

% function [perf, ori, moco, moco_norm, ... 
%     aif_im, aif_moco, aif_cin, aif_cin_Gd, aif_cin_Gd_without_R2Star, aif_cin_Gd_baseline_corrected, ... 
%     aif_cin_all_echo0_signal, aif_cin_all_echo1_signal, aif_cin_all_echo0_signal_after_R2StarCorrection, ...
%     aif_cin_all_echo0_OverPD_after_R2StarCorrection, aif_cin_all_R2Star,  aif_cin_all_R2Star_SLEP, ... 
%     aif_PD, aif_mask, aif_mask_final, aif, aif_baseline_corrected, ... 
%     flow, Ki, PS, Vp, Visf, E, SDMap, Delay, ...
%     BTEX_Flow_all, BTEX_PS_all, BTEX_Vp_all, BTEX_Visf_all, BTEX_cost_all, BTEX_flow_SD_all] = read_in_results(resDir)
% 
%     slc = 0;   
%     for n=1:8
%         
%         filename = ['perf_moco_upsampled_' num2str(n-1) '.hdr'];
%         if(~isFileExist(fullfile(resDir, 'DebugOutput', filename)))
%             break;
%         end        
%         slc = slc + 1;
%     end
%     
%     disp(['Total ' num2str(slc) ' is found ...']);
%     
%     try
%         perf = load_array(resDir, 'CASignal_Perf_', slc);        
%         perf = permute(perf, [1 2 4 3]);
%         perf = flipdim(perf, 2);
%     catch
%         perf = readGTPlusExportImageSeries_Squeeze(resDir, 111);
%     end
% 
%     ori = load_array(resDir, 'perf_', slc);   
%     ori = flipdim(ori, 2);
% 
%     moco = load_array(resDir, 'perf_moco_', slc);   
%     moco = flipdim(moco, 2);
%     
%     moco_norm = load_array(resDir, 'SRNorm_', slc);   
%     moco_norm = flipdim(moco_norm, 2);
% 
%     try
%         aif_im = readGTPlusExportImageSeries_Squeeze(resDir, 1104);
%     catch
%         aif_im = [];
%     end
% 
%     try
%         r1 = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_moco.hdr'));
%         r2 = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_moco_second_echo.hdr'));
%         aif_moco = flipdim(cat(4, r1, r2), 2);
%     catch
%         aif_moco = [];
%     end
% 
%     try
%         aif_cin = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin.hdr'));
%         aif_cin_Gd = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo0_LUTCorrection.hdr'));
%         aif_cin_Gd_without_R2Star = analyze75read(fullfile(resDir, 'DebugOutput', 'cin_all_echo0_without_R2Star_LUTCorrection.hdr'));
%         aif_cin_Gd_baseline_corrected = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_echo0_all_signal_baseline_corrected.hdr'));
%         aif_cin_all_echo0_signal = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo0_signal.hdr'));
%         aif_cin_all_echo1_signal = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo1_signal.hdr'));
%         aif_cin_all_echo0_signal_after_R2StarCorrection = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo0_signal_after_R2StarCorrection.hdr'));
%         aif_cin_all_echo0_OverPD_after_R2StarCorrection = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo0_OverPD_after_R2StarCorrection.hdr'));
%         aif_cin_all_R2Star = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_R2Star.hdr'));
%         aif_cin_all_R2Star_SLEP = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_R2Star_SLEP.hdr'));
%         aif_PD = analyze75read(fullfile(resDir, 'DebugOutput', 'aifPD_for_TwoEcho_T2StartCorrection_0.hdr'));
%         aif_mask = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_LV_mask_for_TwoEcho_T2StartCorrection_0.hdr'));
%         aif_mask_final = analyze75read(fullfile(resDir, 'DebugOutput', 'AifLVMask_after_Picking.hdr'));
% 
%         aif_mask = flipdim(aif_mask, 2);
%         aif_mask_final = flipdim(aif_mask_final, 2);
%         aif_PD = flipdim(aif_PD, 2);
% 
%         aif = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_all_echo0_LUTCorrection.hdr'));
%     catch
%         aif = [];
%     end
% 
%     try
%         aif_baseline_corrected = analyze75read(fullfile(resDir, 'DebugOutput', 'aif_cin_echo0_all_signal_baseline_corrected.hdr'));
%     catch
%         aif_baseline_corrected = [];
%     end
% 
%     try
%         flow = load_array(resDir, 'flow_maps_after_hole_filling_', slc);
%         flow = squeeze(flow(:,:,end,:));
%         flow = flipdim(flow, 2);
%     catch
%         flow = [];
%     end
% 
%     try
%         Ki = load_array(resDir, 'Ki_maps_after_hole_filling_', slc);        
%         Ki = flipdim(Ki, 2);
%     catch
%         Ki = [];
%     end
% 
%     try
%         PS = load_array(resDir, 'PS_maps_after_hole_filling_', slc);
%         PS = squeeze(PS(:,:,end,:));
%         PS = flipdim(PS, 2);
%     catch
%         PS = [];
%     end
% 
%     try
%         Vp = load_array(resDir, 'blood_volume_maps_after_hole_filling_', slc);
%         Vp = squeeze(Vp(:,:,end,:));
%         Vp = flipdim(Vp, 2);
%     catch
%         Vp = [];
%     end
% 
%     try
%         Visf = load_array(resDir, 'interstitial_volume_maps_', slc);
%         Visf = squeeze(Visf(:,:,end,:));
%         Visf = flipdim(Visf, 2);
%     catch
%         Visf = [];
%     end
% 
%     try
%         E = load_array(resDir, 'E_maps_', slc);
%         E = squeeze(E(:,:,end,:));
%         E = flipdim(E, 2);
%     catch
%         E = [];
%     end
% 
%     try
%         for n=1:slc
%             filename = ['BTEX_SD_maps_' num2str(n-1) '_0.hdr'];
%             SDMap(:,:,n) = analyze75read(fullfile(resDir, 'DebugOutput', filename));
%         end
%         SDMap = flipdim(SDMap, 2);
%     catch
%         SDMap = [];
%     end
% 
%     try
%         Delay = load_array(resDir, 'BTEX_res_', slc);
%         Delay = squeeze(Delay(:,:,end,:));
%         Delay = flipdim(Delay, 2);
%     catch
%         Delay = [];
%     end
% 
%     try
%         BTEX_Flow_all = load_array(resDir, 'BTEX_Flow_all_', slc);
%         BTEX_Flow_all = flipdim(BTEX_Flow_all, 2);
%     catch
%         BTEX_Flow_all = [];
%     end
% 
%     try
%         BTEX_PS_all = load_array(resDir, 'BTEX_PS_all_', slc);
%         BTEX_PS_all = flipdim(BTEX_PS_all, 2);
%     catch
%         BTEX_PS_all = [];
%     end
% 
%     try
%         BTEX_Visf_all = load_array(resDir, 'BTEX_Visf_all_', slc);
%         BTEX_Visf_all = flipdim(BTEX_Visf_all, 2);
%     catch
%         BTEX_Visf_all = [];
%     end
% 
%     try
%         BTEX_Vp_all = load_array(resDir, 'BTEX_Vp_all_', slc);
%         BTEX_Vp_all = flipdim(BTEX_Vp_all, 2);
%     catch
%         BTEX_Vp_all = [];
%     end
% 
%     try
%         BTEX_cost_all = load_array(resDir, 'BTEX_cost_all_', slc);
%         BTEX_cost_all = flipdim(BTEX_cost_all, 2);
%     catch
%         BTEX_cost_all = [];
%     end
% 
%     try
%         BTEX_flow_SD_all = load_array(resDir, 'BTEX_flow_SD_all_', slc);
%         BTEX_flow_SD_all = flipdim(BTEX_flow_SD_all, 2);
%     catch
%         BTEX_flow_SD_all = [];
%     end
% end
% 
% function v = load_array(resDir, name, slc)
%     for n=1:slc
%         filename = [name num2str(n-1) '.hdr'];
%         v(:,:,:,n) = analyze75read(fullfile(resDir, 'DebugOutput', filename));
%     end
% end