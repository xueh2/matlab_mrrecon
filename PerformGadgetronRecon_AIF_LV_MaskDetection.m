
function [mask, LV, peak_time] = PerformGadgetronRecon_AIF_LV_MaskDetection(aif_moco, upslope_thres, auc_thres, area_thres, ro_boundary_ratio, e1_boundary_ratio)
% [mask, LV, peak_time] = PerformGadgetronRecon_AIF_LV_MaskDetection(aif_moco)

sampleinterval = 0.5;
sigmas = [0.8 1.2 2.2];
sigmaMeasure = 0.2;
thresGradient = 0.5;

plotFlag = 1;

[upslope, auc, peak_time, foot] = Matlab_gt_perfusion_flow_feature_maps(single(aif_moco));

upslope_sorted = sort(upslope(:));
auc_sorted = sort(auc(:));

num = numel(upslope_sorted);
thres1 = upslope_sorted( floor(upslope_thres*num))

num = numel(auc_sorted);
thres2 = auc_sorted( floor(auc_thres*num))

mask_ind = find(upslope>thres1 & auc>thres2);

mask = upslope;
mask(:) = 0;
mask(mask_ind) = 1;

RO = size(aif_moco,1);
E1 = size(aif_moco,2);
N = size(aif_moco, 3);

sRO = floor(ro_boundary_ratio*RO);
eRO = RO-sRO;

sE1 = floor(e1_boundary_ratio*E1);
eE1 = E1-sE1;

cRO = (sRO+eRO)/2;
cE1 = (sE1+eE1)/2;

a = (eRO-sRO)/2;
b = (eE1-sE1)/2;

for e1=1:E1
    for ro=1:RO

        if(ro<=sRO | ro>=eRO | e1<=sE1 | e1>=eE1)
            mask(ro, e1) = 0;
            continue;
        end

        x = ro - cRO;
        y = e1 - cE1;
        t = (x/a)^2 + (y/b)^2;
        if(t>1)
            mask(ro, e1) = 0;
        end
    end

end

[L, numROI] = bwlabel(mask, 4);        
D=regionprops(L, 'area');

for e1=1:E1
    for ro=1:RO

        if(L(ro,e1)>0)
            if(D(L(ro,e1)).Area<area_thres)
                mask(ro, e1) = 0;
            end
        end

    end
end

opts = statset('Display', 'off');

ind = find(mask>0);

mask_all = repmat(mask, [1 1 N]);
aif_im = aif_moco .* mask_all;
%         figure; imagescn(aif_im, [], [], [], 3);

X = zeros(numel(ind), N);

for n=1:numel(ind)
    [ro, e1] = ind2sub([RO E1], ind(n));
    X(n, :) = aif_im(ro, e1, :);
end

[idx, ctrs] = kmeans(X, 2, 'Distance', 'sqEuclidean', 'Replicates', 10, 'Options', opts);

k_label = mask;
for n=1:numel(ind)
    [ro, e1] = ind2sub([RO E1], ind(n));
    k_label(ro, e1) = idx(n);
end
figure; imagescn(k_label); colormap(jet);

%         closeall

% find the LV curve

baseline = mean(ctrs(:,1:4), 2);
baseline = repmat(baseline, [1 N]);
ctrs = ctrs - baseline;
figure; plot(ctrs'); axis([1 size(ctrs,2) 0 max(ctrs(:))])

[v, ind_c] = max(ctrs');
if(max(ind_c)==ind_c(1))
    LV = ctrs(1,:);
    LV = LV(:);

    LV_ind = 1;
    
    RV = ctrs(2,:);
    RV = RV(:);
else
    LV = ctrs(2,:);
    LV = LV(:);

    LV_ind = 2;
    
    RV = ctrs(1,:);
    RV = RV(:);
end

figure;
hold on
plot(LV);
plot(RV, 'r');
hold off
legend('LV', 'RV');

LV_peak_time = [];
ind_pt = 1;
k_label = mask;
for n=1:numel(ind)
    [ro, e1] = ind2sub([RO E1], ind(n));
    if(idx(n)==LV_ind)
        LV_peak_time(ind_pt) = peak_time(ro, e1);
        ind_pt = ind_pt+1;
    end
end

if(std(LV_peak_time)>3)
    % try again
    
    [idx, ctrs] = kmeans(X, 3, 'Distance', 'sqEuclidean', 'Replicates', 10, 'Options', opts);
    figure; plot(ctrs'); axis([1 size(ctrs,2) 0 max(ctrs(:))])

    k_label = mask;
    for n=1:numel(ind)
        [ro, e1] = ind2sub([RO E1], ind(n));
        k_label(ro, e1) = idx(n);
    end
    figure; imagescn(k_label); colormap(jet);

    baseline = mean(ctrs(:,1:4), 2);
    baseline = repmat(baseline, [1 N]);
    ctrs = ctrs - baseline;
    figure; plot(ctrs'); axis([1 size(ctrs,2) 0 max(ctrs(:))])

    %         closeall

    % find the LV curve
    [v, ind_c] = max(ctrs');   
    [sv, sv_ind] = sort(v);
    
    l2 = sv_ind(end);
    l1 = sv_ind(end-1);
    
    [slope_1, timeToPeak_1, peakTime_1, areaUnderCurve_1, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(ctrs(l1, :)', sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, plotFlag, 'Feature detection of AIF LV signal');
    [slope_2, timeToPeak_2, peakTime_2, areaUnderCurve_2, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(ctrs(l2, :)', sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, plotFlag, 'Feature detection of AIF LV signal');

    if( abs(peakTime_1-peakTime_2)<4)
        l1 = sv_ind(end-2);
    end
    
    ind_c = [ind_c(l1) ind_c(l2)];
    ctrs_2 = [ctrs(l1, :); ctrs(l2, :)];
    
    if(max(ind_c)==ind_c(1))
        LV = ctrs_2(1,:);
        LV = LV(:);

        LV_ind = 1;

        RV = ctrs_2(2,:);
        RV = RV(:);
    else
        LV = ctrs_2(2,:);
        LV = LV(:);

        LV_ind = 2;

        RV = ctrs_2(1,:);
        RV = RV(:);
    end
    
    figure;
    hold on
    plot(LV);
    plot(RV, 'r');
    hold off
    legend('LV', 'RV');

end

[slope, timeToPeak, peakTime, areaUnderCurve, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(RV, sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, plotFlag, 'Feature detection of AIF LV signal');

foot = peakTime - timeToPeak;

foot = floor(foot);

meanLV = mean(LV(1:foot));
LV = LV - meanLV;

[slope_LV, timeToPeak_LV, peakTime_LV, areaUnderCurve_LV, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(LV, sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, plotFlag, 'Feature detection of AIF LV signal');

% then 4 classes
[idx, ctrs] = kmeans(X, 4, 'Distance', 'sqEuclidean', 'Replicates', 10, 'Options', opts);

ctrs = ctrs';

meanCtrs = mean(ctrs(1:foot, :), 1);
meanCtrs = repmat(meanCtrs, [N 1]);

ctrs = ctrs - meanCtrs;
figure; 
for n=1:size(ctrs,2)   
    subplot(size(ctrs,2), 1, n);
    hold on
    plot(LV, 'r--');
    plot(ctrs(:,n)); axis([1 size(ctrs,1) 0 max(ctrs(:))])
    hold off
end

k_label = zeros(size(mask));
for n=1:numel(ind)
    [ro, e1] = ind2sub([RO E1], ind(n));
    k_label(ro, e1) = idx(n);
end
figure; imagescn(k_label); colormap(jet);

[m_ctrs, loc] = max(ctrs)

similarity = [];
for n=1:size(ctrs, 2)
    [r, p] = corrcoef(LV', ctrs(:,n)');
    similarity(n) = r(1,2);

%             similarity(n) = LV' * ctrs(:,n);
end

[sim_sorted, ind_sorted] = sort(abs(similarity))

if(sim_sorted(end)/sim_sorted(end-1)>1.1)
    maxROI = ind_sorted(end);
else
    % use the one with higher peak and bigger
    s1 = ctrs(:, ind_sorted(end));
    num_s1 = numel(find(idx==ind_sorted(end)));
    
    s2 = ctrs(:, ind_sorted(end-1));
    num_s2 = numel(find(idx==ind_sorted(end-1)));

    if(max(s2)>max(s1) & num_s2>=0.8*num_s1)
        maxROI = ind_sorted(end-1);
    else
        maxROI = ind_sorted(end);
    end            
end

k_label = zeros(size(mask));
for n=1:numel(ind)
    [ro, e1] = ind2sub([RO E1], ind(n));
    if(idx(n)==maxROI)
        k_label(ro, e1) = 1;
    end
end
figure; imagescn(k_label);

[L, numROI] = bwlabel(k_label, 4);        
if(numROI>1)
    D=regionprops(L, 'area');

    for e1=1:E1
        for ro=1:RO                
            if(L(ro,e1)>0)
                if(D(L(ro,e1)).Area<max(D.Area))
                    k_label(ro, e1) = 0;
                end
            end                
        end
    end
end

mask = k_label;
peak_time = peakTime_LV;
