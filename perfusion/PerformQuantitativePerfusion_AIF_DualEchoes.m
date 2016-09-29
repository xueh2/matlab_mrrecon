
function flowmaps = PerformQuantitativePerfusion_AIF_DualEchoes(AIF, perf, AIF_mask, AIF_acqTime, perf_mask, perf_acqTime, upsamplingRatio)
% perform the quantitative perfusion computation for dual sequence

%% upsample AIF time stamps to regular grid

% figure; imagescn(AIF, [], [], [], 3);
% figure; imagescn(perf, [], [], [], 3);

interp_method = 'cubic';

REP = size(AIF, 3);
SLC = size(perf, 4);

t_start = AIF_acqTime(1);
t_end = AIF_acqTime(REP);

N = floor(REP*upsamplingRatio);

deltaT = (t_end - t_start) / (N-1);
AIF_time_stamp = t_start:deltaT:t_end;

AIF_upsampled = interp1(AIF_acqTime, permute(AIF, [3 1 2]), AIF_time_stamp, interp_method);
AIF_upsampled = permute(AIF_upsampled, [2 3 1]);

AIF_ori_upsampled = interp1(AIF_acqTime, permute(AIF_ori, [3 1 2]), AIF_time_stamp, interp_method);
AIF_ori_upsampled = permute(AIF_ori_upsampled, [2 3 1]);

% figure; imagescn(AIF_upsampled, [], [], [], 3);
% figure; imagescn(AIF_ori_upsampled, [], [], [], 3);

%% upsample the perf signal
% enforce the border-value boundary condition

RO = size(perf, 1);
E1 = size(perf, 2);

RO_AIF = size(AIF, 1);
E1_AIF = size(AIF, 2);


perf_upsampled = zeros(RO, E1, N, SLC);

for s=1:SLC
       
    perfAcq = perf_acqTime(:,s);
    perfAcq = [AIF_time_stamp(1)-0.01; perfAcq];
    
    perf_slc = zeros(RO, E1, REP+1);
    perf_slc(:,:,1) = perf(:,:,1,s);
    perf_slc(:,:,2:end) = perf(:,:,:,s);
    
    perf_upsampled_slc = interp1(perfAcq, permute(perf_slc, [3 1 2]), AIF_time_stamp, interp_method);
    perf_upsampled_slc = permute(perf_upsampled_slc, [2 3 1]);
    perf_upsampled(:,:,:,s) = perf_upsampled_slc;
end

% figure; imagescn(perf_upsampled, [], [], [], 3);

%% get LV ROI from AIF images

% compute std images
stdAIF = std(AIF_upsampled, 1, 3);
header = CreateFtkHeaderInfo(AIF_upsampled, [1 1 1 1]);
sampleinterval = 1;
sigmas = [3.2 4.0 5.6];
sigmaMeasure = 1.6;
thresGradient = 0.2;
[S, TTP, PT, AUC] = ComputePerfusionParameterMap(AIF_ori_upsampled, header, sampleinterval, sigmas, sigmaMeasure, thresGradient);

ind = find(AIF_mask(:)==0);
if ( ~isempty(ind) )
    S(ind(:)) = 0;
    AUC(ind(:)) = 0;
end

S_sort = sort(S(:));
minGrey = S_sort( floor(0.85*numel(S_sort)) );
maxGrey = S_sort(end);

result_S = zeros(size(S));
result_S(find(S(:)>minGrey)) = 1;

% binNumber = 64;
% relaxationRatio = 0.2;
% maxIter = 3;
% numOfPixelThres = 10;
% [result_S, optimalThres_S] = Otsu_ThresholdRelaxation_UP(S, minGrey, maxGrey, binNumber, relaxationRatio, maxIter, numOfPixelThres);

ind = find(AUC>0);
AUC_sort = sort(AUC(ind(:)));
minGrey = AUC_sort(floor(RO_AIF*E1_AIF*0.05));
maxGrey = AUC_sort(end-floor(RO_AIF*E1_AIF*0.05));
binNumber = 64;
relaxationRatio = 0.2;
maxIter = 3;
numOfPixelThres = 10;
[result_AUC, optimalThres_AUC] = Otsu_ThresholdRelaxation_UP(AUC, minGrey, maxGrey, binNumber, relaxationRatio, maxIter, numOfPixelThres);

mask = result_S.*result_AUC;
mask(find(mask>0)) = 1;
mask(1:floor(0.2*RO_AIF), :) = 0;
mask(floor(0.8*RO_AIF):RO_AIF, :) = 0;
mask(:, 1:floor(0.1*E1_AIF)) = 0;
mask(:, floor(0.9*E1_AIF):E1_AIF) = 0;

PT_masked = PT.*mask;
L = bwlabel(mask, 4);
maxL = max(L(:))

meanPT = -1;
stdmeanPT = -1;
bestL = 0;
for l=1:maxL
    ind = find(L(:)==l);
    if (numel(ind) > 10 )
        PTs = PT(ind(:));
        meanPTCurr = mean(PTs(:))
        stdPTCurr = std(PTs(:))
        if ( stdPTCurr<10 && meanPTCurr > meanPT )
            stdmeanPT = stdPTCurr*meanPTCurr;
            meanPT = meanPTCurr;
            bestL = l;
        end
    end
end

if bestL ~= 0
    mask = zeros(size(S));
    mask(L(:)==bestL) = 1;
end

figure; imagescn(cat(3, AIF_upsampled(:,:,round(meanPT)), mask));
AIF_LV_ROI = mask;

%% get the AIF signal and foot time
ind = find(mask(:)>0);
num = numel(ind(:));

maskAll = repmat(mask, [1 1 N]);
AIF_signal_REP = maskAll.*AIF_upsampled;

AIF_signal = zeros(N, 1);
for r=1:N
   AIF2D = AIF_signal_REP(:,:,r);
   ind2D = find(AIF2D(:)>0);
   v = AIF2D(ind2D(:));
   
   v2 = sort(v);
   nn = numel(v);
   AIF_signal(r) = mean(v2(floor(nn/2):end));
end

[slope, timeToPeak, peakTime, areaUnderCurve, goodFlag] = ...
            PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(AIF_signal, sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, 1, 'AIF signal');

foot = floor(peakTime - timeToPeak);
AIF_foottime = foot;

AIF_signal_baseline_corrected = AIF_signal - mean(AIF_signal(1:foot-1));

perf_upsampled_baseline = mean(perf_upsampled(:,:,1:foot-1,:), 3);

perf_upsampled_baseline = zeros(RO, E1, 1, SLC);
for s=1:SLC  
    for e1=1:E1
        for ro=1:RO
           
            v = 0;
            m = 1;
            for t=1:foot/2
                if (perf_upsampled(ro, e1, t, s)~=-1)
                    v = v+perf_upsampled(ro, e1, t, s);
                    m = m+1;
                end
            end
            
            perf_upsampled_baseline(ro, e1, 1, s) = v/m;
        end
    end    
end

perf_upsampled_baseline_corrected = perf_upsampled - repmat(perf_upsampled_baseline, [1 1 N 1]);

%% perform deconvolution
thres_svd = 1e-4;
orderBSpline = 4;
% numOfInternalControlPoints = 15;
% lambda = 0.0025;

r = zeros(N-foot+1, RO, E1, SLC);
F = zeros(RO, E1, SLC);
F2 = zeros(RO, E1, SLC);
F3 = zeros(RO, E1, SLC);
F4 = zeros(RO, E1, SLC);

cin = AIF_signal_baseline_corrected(foot:end);
for s=1:SLC  
    disp(['----> SLC = ' num2str(s) ' <----']);
    for e1=1:E1
        disp(['      ----> e1 = ' num2str(e1) ' <----']);
        for ro=1:RO
            
%             if (perf_mask(ro, e1)==0 )
%                 continue;
%             end
            
            y_ori = perf(ro, e1, :, s);
            y_ori = squeeze(y_ori);
            
            indY = find(y_ori(:)==-1 );
            if ( ~isempty(indY) | (sum(y_ori(:)) == 0) )
                r(:, ro, e1, s) = -1;
                F(ro, e1, s) = -1;
                continue;
            end
            
            y = perf_upsampled_baseline_corrected(ro, e1, foot:end, s);
            
            y = squeeze(y);
%           imagescn(perf_upsampled_baseline_corrected, [], [], [], 3);
%           figure; hold on; plot(cin, '.--'); plot(y, 'r+'); hold off
%           figure; hold on; plot(cin, '.--'); plot(y, 'r+'); plot(yr,'-x'); plot(yr2, 'g-x'); plot(yr3, 'y-o'); plot(yr4, 'k-*'); hold off; legend('cin', 'y', 'L2 BSpline', 'L1', 'L1 BSpline', 'Fermi');
%           figure; hold on; plot(impulse, '.--'); plot(impulse2, 'r-'); plot(impulse3, 'y-'); plot(impulse4, 'g-'); hold off; legend('L2 BSpline', 'L1', 'L1 BSpline', 'Fermi');

            [impulse, yr] = PerformDeconvolution_Tikhonov_BSpline(cin, y, thres_svd, orderBSpline, numOfInternalControlPoints, lambda);
            [impulse2, yr2] = PerformDeconvolution_L1_Fista(cin, y, lambda);
            [impulse3, yr3] = PerformDeconvolution_L1_Fista_BSpline(cin, y, orderBSpline, numOfInternalControlPoints_L1, lambda);
            
            [v, ind] = PerformDeconvolution_findFirstPeak(impulse3, 1);
            [FF, tau, k, td, yr4, impulse4] = PerformDeconvolution_Fermi(cin, y, v);
            
            r(:, ro, e1, s) = impulse(:);

            [v, ind] = PerformDeconvolution_findFirstPeak(impulse, 0.9);
            F(ro, e1, s) = v;
           
            [v, ind] = PerformDeconvolution_findFirstPeak(impulse2, 0.9);
            F2(ro, e1, s) = v;
                                 
            [v, ind] = PerformDeconvolution_findFirstPeak(impulse3, 0.9);
            F3(ro, e1, s) = v;
            
            % [v, ind] = PerformDeconvolution_findFirstPeak(impulse4, 0.9);
            F4(ro, e1, s) = max(impulse4);
            
        end
    end
end


% convert to ml/min/g
F = F*60*1000/deltaT/1.06;
F2 = F2*60*1000/deltaT/1.06;
F3 = F3*60*1000/deltaT/1.06;
F4 = F4*60*1000/deltaT/1.06;
