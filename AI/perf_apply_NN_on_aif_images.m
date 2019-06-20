
function [lv, rv, probs] = perf_apply_NN_on_aif_images(trainingDataDir, lv_model, lv_rv_model, contourDir)

if(nargin<5)
    contourDir = trainingDataDir;
end
    
lv = [];
rv = [];
probs = [];

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

rE1 = 48;
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
end

lv = analyze75read(fullfile(contourDir, 'LV_Mask'));
probs = analyze75read(fullfile(contourDir, 'prob'));

% ---------------------------------------------------------------------------------------------------------
% lv_rv

if(~isempty(lv_rv_model))
    command = ['gadgetron_CMR_ML_aif_lv_mask --image_file ' fullfile(trainingDataDir, 'aif_for_training') ... 
        ' --prob ' fullfile(contourDir, 'prob') ...
        ' --lv ' fullfile(contourDir, 'LV_Mask') ... 
        ' --rv ' fullfile(contourDir, 'RV_Mask') ... 
        ' --thres 0.5' ...
        ' --rv_thres 0.5' ...
        ' -t lv_rv --model_file ' lv_rv_model ];
    command
    dos(command);    
end

lv = analyze75read(fullfile(contourDir, 'LV_Mask'));
rv = analyze75read(fullfile(contourDir, 'RV_Mask'));
probs = analyze75read(fullfile(contourDir, 'prob'));

