
function PerformGadgetronRecon_AIF_LV_MaskDetection_Cases(perf_cases, resDir, fig_Dir)
% PerformGadgetronRecon_AIF_LV_MaskDetection_Cases(perf_cases, fig_Dir)

startN = 1
endN = size(perf_cases, 1)

upslope_thres = 0.6
upslope_high_thres = 0.9 
auc_thres = 0.9
area_thres = 10
ro_boundary_ratio = 0.15
e1_boundary_ratio = 0.15
max_duration = 40;
sampleinterval = 0.5;
sigmas = [0.8 1.2 2.2];
sigmaMeasure = 0.2;
thresGradient = 0.5;
LV_picked = 0.15;
plotFlag = 1;

failed_cases = [];

for ii=startN:endN
    ii
    try
        disp([num2str(ii) ' of ' num2str(endN) ' ... ']);
        [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_stress] = parseSavedISMRMRD(perf_cases{ii,2});

        % stress
        stress_dir = fullfile(resDir, study_dates, perf_cases{ii,2})

        aif_moco_stress = analyze75read(fullfile(stress_dir, 'DebugOutput', 'aif_upsampled.hdr'));
        aif_pd_stress = analyze75read(fullfile(stress_dir, 'DebugOutput', 'aifPD_for_TwoEcho_T2StartCorrection_0.hdr'));
        aif_moco_stress_mask = analyze75read(fullfile(stress_dir, 'DebugOutput', 'aif_LV_mask.hdr'));

        [mask_stress, LV, peakTime_LV_stress, mask_RV_stress] = PerformGadgetronRecon_AIF_LV_MaskDetection(aif_moco_stress, upslope_thres, auc_thres, area_thres, max_duration, ro_boundary_ratio, e1_boundary_ratio);
       
        ind = find(mask_stress>0);
        ind_RV = find(mask_RV_stress>0);
        
        aif_slc = aif_moco_stress(:,:, floor(peakTime_LV_stress));
        
        aif_slc_mask = aif_slc;
        aif_slc_mask(ind) = max(aif_slc(:))+1;
        
        [L, numroi] = bwlabel(mask_stress, 4);
        C = regionprops(L, 'Centroid');
        R = regionprops(L, 'EquivDiameter');
        
        RO = size(aif_slc_mask, 1);
        E1 = size(aif_slc_mask, 2);
        N = size(aif_moco_stress, 3);
        
        ind = find(aif_moco_stress_mask>0);
        aif_slc_mask_old = aif_slc;
        aif_slc_mask_old(ind) = max(aif_slc(:))+1;
        
        h = figure('Name',perf_cases{ii,3},'NumberTitle','off'); imagescn(cat(3, aif_slc, aif_slc_mask, aif_slc), [], [1 3], [10]);
        posh = get(h, 'Position');
        set(h, 'Position', [posh(1) 12 posh(3) posh(4)]);
        hold on
        z = R.EquivDiameter/2*exp((0:.1:2*pi)*i);
        plot(real(z)+C.Centroid(1), imag(z)+C.Centroid(2), 'r-', 'LineWidth', 2);
        hold off
        saveas(h, fullfile(fig_Dir, perf_cases{ii,2}), 'fig');
        ff = getframe(h);
        [X, map] = frame2im(ff);
        imwrite(X, fullfile(fig_Dir, [perf_cases{ii,2} '.bmp']), 'BMP');
        
        % compute LV signals
        [slope_LV, timeToPeak_LV, peakTime_LV, areaUnderCurve_LV, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(LV, sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, plotFlag, 'Feature detection of AIF LV signal');
        
        numLV = numel(ind);
        
        X = zeros(N, numLV);
        for nn=1:numLV
            for f=1:N
                af = aif_moco_stress(:,:,f);
                a = af(ind(:));
                X(f, :) = a(:);
            end
        end
        
        peakX = zeros(numLV, 1);
        for nn=1:numLV
            peakX(nn,1) = X(floor(peakTime_LV), nn) + X(ceil(peakTime_LV), nn);
        end
        
        [peakX_sorted, ind_peakX] = sort(peakX, 1, 'descend');
        
        selected_N = floor(numLV*LV_picked);
        if(selected_N<1) 
            selected_N = 1;
        end
        
        X_selected = X(:, ind_peakX(1:selected_N));
        
        X_mean_all = mean(X, 2);
        X_mean = mean(X_selected, 2);
        
        [V,D] = KL_Eigen(reshape(X_selected, [N 1 selected_N]));
        
        X_eigen = X_selected*V;
        X_eigen = X_eigen / sqrt(selected_N);
        
        X_eigen_mean = X_eigen(:, end);
        
        if( X_mean(floor(peakTime_LV))*X_eigen_mean(floor(peakTime_LV))<0 )
            X_eigen_mean = X_eigen_mean * -1;
        end
        
        % compute RV signal       
        RV = zeros(N, 1);
        for n=1:N
            af = aif_moco_stress(:,:,n);
            RV(n) = mean(af(ind_RV(:)));
        end
        
        % plot
        
        t = 0.5 * [0:N-1];
        h2 = figure('Name',perf_cases{ii,2},'NumberTitle','off', 'Position', [40 300 1500 700]);
        hold on
        plot(t, X_eigen_mean, 'r', 'LineWidth', 4);
        plot(t, X_mean, 'b', 'LineWidth', 2);
        plot(t, X_mean_all, 'k', 'LineWidth', 2);
        plot(t, RV, 'b--', 'LineWidth', 3);
    	hold off
        legend('Eigen mean with top 15%', 'Norm mean with top 15%', 'Norm mean with all', 'RV mean with all', 'Location', 'NorthWest');
        xlabel('seconds')
        ylabel('AIF intensity')
        box on
        grid on
        
        saveas(h2, fullfile(fig_Dir, [perf_cases{ii,2} '_eigenMean']), 'fig');
        ff = getframe(h2);
        [X, map] = frame2im(ff);
        imwrite(X, fullfile(fig_Dir, [perf_cases{ii,2} '_eigenMean.bmp']), 'BMP');       
        
        % rest
        if(size(perf_cases, 2)>=3)
            
            closeall
            clear C R L mask_stress z ff X map h mask LV peakTime_LV aif_slc_mask X_eigen X RV
            
            [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_rest] = parseSavedISMRMRD(perf_cases{ii, 3});
            rest_dir = fullfile(resDir, study_dates, perf_cases{ii,3})

            aif_moco_rest = analyze75read(fullfile(rest_dir, 'DebugOutput', 'aif_upsampled.hdr'));
            aif_pd_rest = analyze75read(fullfile(rest_dir, 'DebugOutput', 'aifPD_for_TwoEcho_T2StartCorrection_0.hdr'));
            aif_moco_rest_mask = analyze75read(fullfile(rest_dir, 'DebugOutput', 'aif_LV_mask.hdr'));

            [mask, LV, peakTime_LV, mask_RV_rest] = PerformGadgetronRecon_AIF_LV_MaskDetection(aif_moco_rest, upslope_thres, auc_thres, area_thres, max_duration, ro_boundary_ratio, e1_boundary_ratio);

            ind = find(mask>0);
            ind_RV = find(mask_RV_rest>0);
            
            aif_slc = aif_moco_rest(:,:, floor(peakTime_LV));

            aif_slc_mask = aif_slc;
            aif_slc_mask(ind) = max(aif_slc(:))+1;

            [L, numroi] = bwlabel(mask, 4);
            C = regionprops(L, 'Centroid');
            R = regionprops(L, 'EquivDiameter');

            RO = size(aif_slc_mask, 1);
            E1 = size(aif_slc_mask, 2);
        
            ind = find(aif_moco_rest_mask>0);
            aif_slc_mask_old = aif_slc;
            aif_slc_mask_old(ind) = max(aif_slc(:))+1;

            h = figure('Name',perf_cases{ii,3},'NumberTitle','off'); imagescn(cat(3, aif_slc, aif_slc_mask, aif_slc), [], [1 3], [10]);
            posh = get(h, 'Position');
            set(h, 'Position', [posh(1) 12 posh(3) posh(4)]);
            hold on
            z = R.EquivDiameter/2*exp((0:.1:2*pi)*i);
            plot(real(z)+C.Centroid(1), imag(z)+C.Centroid(2), 'r-', 'LineWidth', 2);
            hold off
            
            saveas(h, fullfile(fig_Dir, perf_cases{ii,3}), 'fig');
            ff = getframe(h);
            [X, map] = frame2im(ff);
            imwrite(X, fullfile(fig_Dir, [perf_cases{ii,3} '.bmp']), 'BMP');
            
            % eigen mean
            % compute LV signals
            [slope_LV, timeToPeak_LV, peakTime_LV, areaUnderCurve_LV, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(LV, sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, plotFlag, 'Feature detection of AIF LV signal');

            numLV = numel(ind);

            X = zeros(N, numLV);
            for nn=1:numLV
                for f=1:N
                    af = aif_moco_rest(:,:,f);
                    a = af(ind(:));
                    X(f, :) = a(:);
                end
            end

            peakX = zeros(numLV, 1);
            for nn=1:numLV
                peakX(nn,1) = X(floor(peakTime_LV), nn) + X(ceil(peakTime_LV), nn);
            end

            [peakX_sorted, ind_peakX] = sort(peakX, 1, 'descend');

            selected_N = floor(numLV*LV_picked);
            if(selected_N<1) 
                selected_N = 1;
            end

            X_selected = X(:, ind_peakX(1:selected_N));

            X_mean_all = mean(X, 2);
            X_mean = mean(X_selected, 2);

            [V,D] = KL_Eigen(reshape(X_selected, [N 1 selected_N]));

            X_eigen = X_selected*V;
            X_eigen = X_eigen / sqrt(selected_N);

            X_eigen_mean = X_eigen(:, end);

            if( X_mean(floor(peakTime_LV))*X_eigen_mean(floor(peakTime_LV))<0 )
                X_eigen_mean = X_eigen_mean * -1;
            end

            % compute RV signal       
            RV = zeros(N, 1);
            for n=1:N
                af = aif_moco_rest(:,:,n);
                RV(n) = mean(af(ind_RV(:)));
            end

            % plot

            t = 0.5 * [0:N-1];
            h2 = figure('Name',perf_cases{ii,3},'NumberTitle','off', 'Position', [40 300 1500 700])
            hold on
            plot(t, X_eigen_mean, 'r', 'LineWidth', 4);
            plot(t, X_mean, 'b', 'LineWidth', 2);
            plot(t, X_mean_all, 'k', 'LineWidth', 2);
            plot(t, RV, 'b--', 'LineWidth', 3);
            hold off
            legend('Eigen mean with top 15%', 'Norm mean with top 15%', 'Norm mean with all', 'RV mean with all', 'Location', 'NorthWest');
            xlabel('seconds')
            ylabel('AIF intensity')
            box on
            grid on

            saveas(h2, fullfile(fig_Dir, [perf_cases{ii,3} '_eigenMean']), 'fig');
            ff = getframe(h2);
            [X, map] = frame2im(ff);
            imwrite(X, fullfile(fig_Dir, [perf_cases{ii,3} '_eigenMean.bmp']), 'BMP');
        end
        
%         pause
        clear C R L mask_stress z ff X map h mask LV peakTime_LV aif_slc_mask X_eigen X_mean RV
        closeall
    catch
    end
end
