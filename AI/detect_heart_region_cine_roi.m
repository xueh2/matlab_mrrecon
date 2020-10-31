function [mask, ro_s, ro_e, e1_s, e1_e, h_mask] = detect_heart_region_cine_roi(data, thres, NN_RO, NN_E1, use_moco, regularization_hilbert_strength)
% detect heart region
% [mask, s_ro, e_ro, s_e1, e_e1] = detect_heart_region_cine_roi(data, thres, NN_RO, NN_E1)

if(nargin<5)
    use_moco = 0;
end

if(nargin<6)
    regularization_hilbert_strength = R0/8;
end

RO = size(data,1);
E1 = size(data,2);
PHS = size(data,3);
SLC = size(data, 4);

strategy = 'FixedReference'; 
dissimilarity = 'LocalCCR'; 
level = 4; 
max_iter_num_pyramid_level = [32 64 100 100]; 
LocalCCR_sigmaArg = 2.0; 
BidirectionalReg = 1; 
dissimilarity_thres = 1e-6; 
div_num = 3; 
inverse_deform_enforce_iter = 10; 
inverse_deform_enforce_weight = 0.5; 
DivergenceFreeReg = 1; 
DebugFolder = []; 
verbose = 0; 
keyframe = round(PHS/2);
if(use_moco)
    [dx, dy, warped, dxInv, dyInv] = Matlab_gt_deformation_field_reg_2D_series(data, keyframe, strategy, dissimilarity, level, max_iter_num_pyramid_level, regularization_hilbert_strength, LocalCCR_sigmaArg, BidirectionalReg, dissimilarity_thres, div_num, inverse_deform_enforce_iter, inverse_deform_enforce_weight, DivergenceFreeReg, DebugFolder, verbose); 
    data = warped;
end

p = std(data, [], 3);
pp = squeeze(sum(p, 4));

sRO = round(0.2*RO);
eRO = round(0.8*RO);
sE1 = round(0.15*E1);
eE1 = round(0.85*E1);

pp(1:sRO, :) = 0;
pp(eRO:RO, :) = 0;
pp(:, 1:sE1) = 0;
pp(:, eE1:E1) = 0;

%figure; imagescn(pp)
[N, X] = hist(pp(:), 100);
ind = cumsum(N)/sum(N);
X_ind = find(ind>thres);
max_p = X(X_ind(1));
ind = find(pp>max_p);
pp2 = zeros(size(pp));
pp2(ind) = 1;
%figure; imagescn(pp2)

if(1)
    L = bwlabel(pp2,8);
    s  = regionprops(L, 'Area');

    max_label = 1;
    max_area = 0;
    valid_labels = [];
    for ii=1:size(s, 1)
        if(s(ii).Area>max_area)
            max_area = s(ii).Area;
            max_label = ii;
        end
    end
    for ii=1:size(s, 1)
        if(s(ii).Area>0.75*max_area)
            valid_labels = [valid_labels; ii];
        end
    end

    ind = find(L==max_label);
    pp2 = zeros(RO, E1);
    pp2(ind) = 1;
    
    for ii=1:numel(valid_labels)
        ind = find(L==valid_labels(ii));
        pp2(ind) = 1;
    end
end

mask = pp2;

[X, Y] = find(pp2==1);

c_ro = round(mean(X));
c_e1 = round(mean(Y));

if (isnan(c_ro))
    c_ro = RO/2;
end
if (isnan(c_e1))
    c_e1 = E1/2;
end

h_mask = -1;
h_mask = figure('Visible', 'off'); imagescn(pp2)
hold on
plot(c_e1, c_ro, '+');
hold off

%figure; 
sp = squeeze(sum(p, 4));
sp(find(pp2==1)) = max(sp(:))+1;
%imagescn(sp);

ro_s = c_ro - NN_RO/2;
if(ro_s<1)
    ro_s = 1;
end
ro_e = ro_s + NN_RO -1;
if(ro_e>RO)
    ro_e = RO;
    ro_s = ro_e - NN_RO + 1;
end

e1_s = c_e1 - NN_E1/2;
if(e1_s<1)
    e1_s = 1;
end
e1_e = e1_s + NN_E1 -1;
if(e1_e>E1)
    e1_e = E1;
    e1_s = e1_e - NN_E1 + 1;
end
