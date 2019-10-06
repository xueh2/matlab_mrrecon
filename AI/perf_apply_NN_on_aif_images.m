
function [lv, rv, probs, lv_3C, rv_3C, probs_3C] = perf_apply_NN_on_aif_images(trainingDataDir, lv_model, lv_rv_model, rE1, contourDir)

if(nargin<6)
    contourDir = trainingDataDir;
end
    
lv = [];
rv = [];
probs = [];

lv_3C = [];
rv_3C = [];
probs_3C = [];

try
    tic
    aif = readNPY(fullfile(trainingDataDir, 'aif_scc.npy'));
    size(aif)
    toc
catch
    tic
    aif = readNPY(fullfile(trainingDataDir, 'aif.npy'));
    size(aif)
    toc
end

RO = size(aif,1);
E1 = size(aif,2);
N = size(aif,3);

oriRO = RO;
oriE1 = E1;

% if(size(aif, 1)==128)
%     aif = Matlab_gt_resize_2D_image(double(aif), RO/2, E1/2, 3);
%     RO = size(aif,1);
%     E1 = size(aif,2);
%     N = size(aif,3);
% end

permute_im = 0;
if(size(aif, 2)==128 | size(aif, 2)==64)
    aif = permute(aif, [2, 1, 3]);
    RO = size(aif,1);
    E1 = size(aif,2);
    N = size(aif,3);
    permute_im = 1;
end

% rE1 = 48;
rN = 64;

if(N>rN)
    aif = aif(:,:,1:rN);
end

if(N<rN)
    aif2 = zeros(RO, E1, rN);
    aif2(:,:,1:N) = aif;
    aif2(:,:,N+1:end) = repmat(aif(:,:,end), [1 1 rN-N]);
    
    aif = aif2;
end

if(E1>rE1)
    oE1 = (E1-rE1)/2;
    aif = aif(:, oE1:oE1+rE1-1, :);
end
if(rE1>E1)
    aif2 = zeros(RO, rE1, rN);
    oE1 = (rE1-E1)/2;
    aif2(:,oE1:oE1+E1-1,:) = aif;
    aif = aif2;
end

header = CreateGtImageHeader(aif);
Matlab_gt_write_analyze(single(aif), header, fullfile(trainingDataDir, 'aif_for_training'));
        
% ---------------------------------------------------------------------------------------------------------
% lv 

if(~isempty(lv_model))
    command = ['gadgetron_CMR_ML_aif_lv_mask --image_file ' fullfile(trainingDataDir, 'aif_for_training') ... 
        ' --prob ' fullfile(contourDir, 'prob') ...
        ' --lv ' fullfile(contourDir, 'LV_Mask') ... 
        ' --thres 0.5' ...
        ' -t lv --model_file ' lv_model ];
    command
    dos(command); 
    
    lv = analyze75read(fullfile(contourDir, 'LV_Mask'));
    probs = analyze75read(fullfile(contourDir, 'prob'));
    
    if(size(lv, 2)<E1)
        lv = zpad(lv, RO, E1);
    end
    
    if(oriRO>RO)
        lv = Matlab_gt_resize_2D_image(double(lv), oriRO, oriE1, 3);
        probs = Matlab_gt_resize_2D_image(double(probs), oriRO, oriE1, 3);
        
        ind = find(lv>0.5);
        lv(:) = 0;
        lv(ind) = 1;
    end
end

% ---------------------------------------------------------------------------------------------------------
% lv_rv

if(~isempty(lv_rv_model))
    command = ['gadgetron_CMR_ML_aif_lv_mask --image_file ' fullfile(trainingDataDir, 'aif_for_training') ... 
        ' --prob ' fullfile(contourDir, 'prob_3C') ...
        ' --lv ' fullfile(contourDir, 'LV_Mask_3C') ... 
        ' --rv ' fullfile(contourDir, 'RV_Mask_3C') ... 
        ' --thres 0.5' ...
        ' --rv_thres 0.5' ...
        ' -t lv_rv --model_file ' lv_rv_model ];
    command
    dos(command);    
    
    lv_3C = analyze75read(fullfile(contourDir, 'LV_Mask_3C'));
    rv_3C = analyze75read(fullfile(contourDir, 'RV_Mask_3C'));
    probs_3C = analyze75read(fullfile(contourDir, 'prob_3C'));
    
    if(size(lv, 2)<E1)
        lv_3C = zpad(lv_3C, RO, E1);
        rv_3C = zpad(rv_3C, RO, E1);
    end
    
    if(oriRO>RO)
        lv_3C = Matlab_gt_resize_2D_image(double(lv_3C), oriRO, oriE1, 3);
        rv_3C = Matlab_gt_resize_2D_image(double(rv_3C), oriRO, oriE1, 3);
        probs_3C = Matlab_gt_resize_2D_image(double(probs_3C), oriRO, oriE1, 3);
        
        ind = find(lv_3C>0.5);
        lv_3C(:) = 0;
        lv_3C(ind) = 1;
    end
end
