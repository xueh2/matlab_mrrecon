
function [mask, LV, peak_time, mask_RV] = PerformGadgetronRecon_AIF_LV_MaskDetection(aif_moco, upslope_thres, auc_thres, area_thres, max_duration, ro_boundary_ratio, e1_boundary_ratio)
% [mask, LV, peak_time, mask_RV] = PerformGadgetronRecon_AIF_LV_MaskDetection(aif_moco, upslope_thres, auc_thres, area_thres, max_duration, ro_boundary_ratio, e1_boundary_ratio)

sampleinterval = 0.5;
sigmas = [0.8 1.2 2.2];
sigmaMeasure = 0.2;
thresGradient = 0.5;

RO = size(aif_moco,1);
E1 = size(aif_moco,2);
N = size(aif_moco, 3);

figure; imagescn(aif_moco, [], [], [], 3);

if(N>max_duration/0.5)
    N = floor(max_duration/0.5);
    aif_moco = aif_moco(:,:,1:N);
end

mask = ones(RO, E1);

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

ind_mask = find(mask>0);

plotFlag = 1;

[upslope, auc, peak_time, foot] = Matlab_gt_perfusion_flow_feature_maps(single(aif_moco));


upslope_sorted = sort(upslope(ind_mask));
auc_sorted = sort(auc(ind_mask));

num = numel(upslope_sorted);
thres1 = upslope_sorted( floor(upslope_thres*num))

num = numel(auc_sorted);
thres2 = auc_sorted( floor(auc_thres*num))

mask_ind = find(upslope>thres1 & auc>thres2 & mask>0);

mask = upslope;
mask(:) = 0;
mask(mask_ind) = 1;

figure; imagescn(mask);

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

figure; imagescn(mask);

[L_2nd, numROI_2nd] = bwlabel(mask, 4);        
D_2nd=regionprops(L_2nd, 'area');
C_2nd=regionprops(L_2nd, 'centroid');
figure; imagescn(L_2nd); colormap(jet);

if(numROI_2nd==2)
    % if one ROI is much smaller than the other, mask it off
    SA(1) = D_2nd(1).Area;
    SA(2) = D_2nd(2).Area;
    
    ind_1 = find(L_2nd==1);
    ind_2 = find(L_2nd==2);
    
    d1 = C_2nd(1).Centroid(:) - [E1/2 RO/2]';
    d1 = norm(d1);
    
    d2 = C_2nd(2).Centroid(:) - [E1/2 RO/2]';
    d2 = norm(d2);
    
    if(SA(2)/SA(1)>2.5 & d2<d1)
        % mask off 1
        ind_2nd = find(L_2nd==1);
        mask(ind_2nd) =  0;
    elseif(SA(1)/SA(2)>2.5 & d1<d2)
        % mask off 2
        ind_2nd = find(L_2nd==2);
        mask(ind_2nd) =  0;        
    end
end

if(numROI_2nd==3)
    
    SA(1) = D_2nd(1).Area;
    SA(2) = D_2nd(2).Area;
    SA(3) = D_2nd(3).Area;
    
    ind_1 = find(L_2nd==1);
    ind_2 = find(L_2nd==2);
    ind_3 = find(L_2nd==3);
    
    d1 = C_2nd(1).Centroid(:) - [E1/2 RO/2]';
    d1 = norm(d1);
    
    d2 = C_2nd(2).Centroid(:) - [E1/2 RO/2]';
    d2 = norm(d2);
    
    d3 = C_2nd(3).Centroid(:) - [E1/2 RO/2]';
    d3 = norm(d3);
    
    [md, vd] = max([d1 d2 d3]);
    
    if(md>0.4*norm([RO/2 E1/2]))
        mask_off_ind = find(L_2nd==vd);
        mask(mask_off_ind) = 0;
        figure; imagescn(mask); title('mask after centroid masking');
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

[slope_1, timeToPeak_1, peakTime_1, areaUnderCurve_1, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(LV, sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, plotFlag, 'Feature detection of AIF LV signal');
[slope_2, timeToPeak_2, peakTime_2, areaUnderCurve_2, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(RV, sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, plotFlag, 'Feature detection of AIF LV signal');

foot_LV = peakTime_1 - timeToPeak_1;
foot_RV = peakTime_2 - timeToPeak_2;

if(1.25*foot_LV < foot_RV & max(LV)>1.1*max(RV))
    % swap LV and RV signal
    temp_s = LV;
    LV = RV;
    RV = temp_s;
    
    if(LV_ind==1)
        LV_ind = 2;
    else
        LV_ind = 1;
    end
    
    figure;
    hold on
    plot(LV);
    plot(RV, 'r');
    hold off
    legend('LV', 'RV');
    title('LV-RV, after swap');
end

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

if(std(LV_peak_time)>3 | (abs(peakTime_1-peakTime_2)<2))
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

    figure; plot(ctrs'); axis([1 size(ctrs,2) 0 max(ctrs(:))]); legend('1', '2', '3');

    %         closeall

    % find the LV curve
    [v, ind_c] = max(ctrs');   
    [sv, sv_ind] = sort(v);
    
    l2 = sv_ind(end);
    l1 = sv_ind(end-1);
    l3 = sv_ind(end-2);
    
    [slope_1, timeToPeak_1, peakTime_1, areaUnderCurve_1, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(ctrs(l1, :)', sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, plotFlag, 'Feature detection of AIF LV signal');
    [slope_2, timeToPeak_2, peakTime_2, areaUnderCurve_2, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(ctrs(l2, :)', sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, plotFlag, 'Feature detection of AIF LV signal');
    [slope_3, timeToPeak_3, peakTime_3, areaUnderCurve_3, goodFlag] = PerfusionParameterEstimation_GaussianSmoothing_RegionMerge(ctrs(l3, :)', sampleinterval, sigmas, sigmaMeasure, 1, thresGradient, plotFlag, 'Feature detection of AIF LV signal');

    foot1 = peakTime_1 - timeToPeak_1;
    foot2 = peakTime_2 - timeToPeak_2;
    
    peakTime_l1 = peakTime_1;
    peakTime_l2 = peakTime_2;
    
    std_baseline1 = std(ctrs(l1,2:floor(foot1-2)), [], 2);
    std_baseline2 = std(ctrs(l2,2:floor(foot2-2)), [], 2);
        
    if( abs(peakTime_1-peakTime_2)<6 | abs(ind_c(l2)-ind_c(l1))<6 | std_baseline1>5*std_baseline2)
        l1 = sv_ind(end-2);
        peakTime_l1 = peakTime_3;
    end
    
    if( abs(peakTime_1-peakTime_2)<4 & abs(peakTime_1-peakTime_3)<4 & abs(peakTime_2-peakTime_3)<4)
        % three curves are all LV
        ind_c = [ind_c(l2) ind_c(l2)];
        ctrs_2 = [ctrs(l2, :); ctrs(l2, :)];
        peak_time_c = [peakTime_l2 peakTime_l2];
    else    
        ind_c = [ind_c(l1) ind_c(l2)];
        ctrs_2 = [ctrs(l1, :); ctrs(l2, :)];
        peak_time_c = [peakTime_l1 peakTime_l2];
    end
    
    %if(max(ind_c)==ind_c(1))
    if(max(peak_time_c)==peak_time_c(1))
        
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
[max_ctrs, ind_RV] = max(m_ctrs);

similarity = [];
for n=1:size(ctrs, 2)
    [r, p] = corrcoef(LV', ctrs(:,n)');
    similarity(n) = r(1,2);

%             similarity(n) = LV' * ctrs(:,n);
end

% ind_one = find(similarity==1);
% if(~isempty(ind_one))
%     similarity(ind_one) = 0;
% end
% 
[sim_sorted, ind_sorted] = sort(abs(similarity))

if(sim_sorted(end)/sim_sorted(end-1)>1.1)
    maxROI = ind_sorted(end);
else
    % use the one with higher peak and bigger
    s1 = ctrs(:, ind_sorted(end));
    num_s1 = numel(find(idx==ind_sorted(end)));
    max_s1 = max(s1);
    
    s2 = ctrs(:, ind_sorted(end-1));
    num_s2 = numel(find(idx==ind_sorted(end-1)));
    max_s2 = max(s2);
    
    if(max_s1/max_s2>1.4)
        maxROI = ind_sorted(end);
    elseif(max_s2/max_s1>1.4)
        maxROI = ind_sorted(end-1);
    else
        if(max(s2)>max(s1) & num_s2>=0.5*num_s1)
            maxROI = ind_sorted(end-1);
        else
            maxROI = ind_sorted(end);
        end            
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
    for rr=1:numROI
        max_D(rr) = D(rr).Area;
    end
    max_D = max(max_D);
    
    for e1=1:E1
        for ro=1:RO                
            if(L(ro,e1)>0)
                if(D(L(ro,e1)).Area<max_D)
                    k_label(ro, e1) = 0;
                end
            end                
        end
    end
end

mask = k_label;
peak_time = peakTime_LV;

% RV
mask_RV = zeros(size(mask));
for n=1:numel(ind)
    [ro, e1] = ind2sub([RO E1], ind(n));
    if(idx(n)==ind_RV)
        mask_RV(ro, e1) = 1;
    end
end
figure; imagescn(mask_RV);
