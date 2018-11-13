
function [endo, epi, endo_epi_rv_rvi] = perf_apply_NN_on_images(trainingDataDir, endo_model, epi_model, endo_epi_rv_rvi_model, contourDir)

if(nargin<5)
    contourDir = trainingDataDir;
end
    
endo = [];
epi = [];
endo_epi_rv_rvi = [];

tic
Gd_resized_training = readNPY(fullfile(trainingDataDir, 'Gd_resized_training.npy'));
size(Gd_resized_training)
toc

if(~isFileExist(fullfile(trainingDataDir, 'Gd_resized_training_1.img')))
    SLC = size(Gd_resized_training, 4);
    for slc=1:SLC
        header = CreateGtImageHeader(Gd_resized_training(:,:,:,slc));
        Matlab_gt_write_analyze(single(Gd_resized_training(:,:,:,slc)), header, fullfile(trainingDataDir, ['Gd_resized_training_' num2str(slc)]));
    end
end

% ---------------------------------------------------------------------------------------------------------
% endo 

if(~isempty(endo_model))
    command = ['gadgetron_CMR_ML_perf_seg_sax -i ' fullfile(trainingDataDir, 'Gd_resized_training_1') ... 
        ' -s ' fullfile(contourDir, 'seg_endo_0') ...
        ' -p ' fullfile(contourDir, 'prob_endo_0') ...
        ' --endo ' fullfile(contourDir, 'endo_0') ... 
        ' -t endo -f ' endo_model ];
    command
    dos(command);

    command = ['gadgetron_CMR_ML_perf_seg_sax -i ' fullfile(trainingDataDir, 'Gd_resized_training_2') ... 
        ' -s ' fullfile(contourDir, 'seg_endo_1') ...
        ' -p ' fullfile(contourDir, 'prob_endo_1') ...
        ' --endo ' fullfile(contourDir, 'endo_1') ... 
        ' -t endo -f ' endo_model ];
    command
    dos(command);

    command = ['gadgetron_CMR_ML_perf_seg_sax -i ' fullfile(trainingDataDir, 'Gd_resized_training_3') ... 
        ' -s ' fullfile(contourDir, 'seg_endo_2') ...
        ' -p ' fullfile(contourDir, 'prob_endo_2') ...
        ' --endo ' fullfile(contourDir, 'endo_2') ... 
        ' -t endo -f ' endo_model ];
    command
    dos(command);
end

endo.seg0 = analyze75read(fullfile(contourDir, 'seg_endo_0'));
endo.seg1 = analyze75read(fullfile(contourDir, 'seg_endo_1'));
endo.seg2 = analyze75read(fullfile(contourDir, 'seg_endo_2'));

endo.prob0 = analyze75read(fullfile(contourDir, 'prob_endo_0'));
endo.prob1 = analyze75read(fullfile(contourDir, 'prob_endo_1'));
endo.prob2 = analyze75read(fullfile(contourDir, 'prob_endo_2'));

endo.C0 = analyze75read(fullfile(contourDir, 'endo_0'));
endo.C1 = analyze75read(fullfile(contourDir, 'endo_1'));
endo.C2 = analyze75read(fullfile(contourDir, 'endo_2'));              


% ---------------------------------------------------------------------------------------------------------
% epi 
if(~isempty(epi_model))
    command = ['gadgetron_CMR_ML_perf_seg_sax -i ' fullfile(trainingDataDir, 'Gd_resized_training_1') ... 
        ' -s ' fullfile(contourDir, 'seg_epi_0') ...
        ' -p ' fullfile(contourDir, 'prob_epi_0') ...
        ' --epi ' fullfile(contourDir, 'epi_0') ... 
        ' -t epi -f ' epi_model ];
    command
    dos(command);

    command = ['gadgetron_CMR_ML_perf_seg_sax -i ' fullfile(trainingDataDir, 'Gd_resized_training_2') ... 
        ' -s ' fullfile(contourDir, 'seg_epi_1') ...
        ' -p ' fullfile(contourDir, 'prob_epi_1') ...
        ' --epi ' fullfile(contourDir, 'epi_1') ... 
        ' -t epi -f ' epi_model ];
    command
    dos(command);

    command = ['gadgetron_CMR_ML_perf_seg_sax -i ' fullfile(trainingDataDir, 'Gd_resized_training_3') ... 
        ' -s ' fullfile(contourDir, 'seg_epi_2') ...
        ' -p ' fullfile(contourDir, 'prob_epi_2') ...
        ' --epi ' fullfile(contourDir, 'epi_2') ... 
        ' -t epi -f ' epi_model ];
    command
    dos(command);
end    
epi.seg0 = analyze75read(fullfile(contourDir, 'seg_epi_0'));
epi.seg1 = analyze75read(fullfile(contourDir, 'seg_epi_1'));
epi.seg2 = analyze75read(fullfile(contourDir, 'seg_epi_2'));

epi.prob0 = analyze75read(fullfile(contourDir, 'prob_epi_0'));
epi.prob1 = analyze75read(fullfile(contourDir, 'prob_epi_1'));
epi.prob2 = analyze75read(fullfile(contourDir, 'prob_epi_2'));

epi.C0 = analyze75read(fullfile(contourDir, 'epi_0'));
epi.C1 = analyze75read(fullfile(contourDir, 'epi_1'));
epi.C2 = analyze75read(fullfile(contourDir, 'epi_2'));


% ---------------------------------------------------------------------------------------------------------
% endo_epi_rv_rvi 
if(~isempty(endo_epi_rv_rvi_model))
    command = ['gadgetron_CMR_ML_perf_seg_sax -i ' fullfile(trainingDataDir, 'Gd_resized_training_1') ... 
        ' -s ' fullfile(contourDir, 'seg_endo_epi_rv_rvi_0') ...
        ' -p ' fullfile(contourDir, 'prob_endo_epi_rv_rvi_0') ...
        ' --endo ' fullfile(contourDir, 'endo_epi_rv_rvi_endo_0') ... 
        ' --epi ' fullfile(contourDir, 'endo_epi_rv_rvi_epi_0') ... 
        ' --rv ' fullfile(contourDir, 'endo_epi_rv_rvi_rv_0') ... 
        ' -t endo_epi_rv_rvi -f ' endo_epi_rv_rvi_model ];
    command
    dos(command);

    command = ['gadgetron_CMR_ML_perf_seg_sax -i ' fullfile(trainingDataDir, 'Gd_resized_training_2') ... 
        ' -s ' fullfile(contourDir, 'seg_endo_epi_rv_rvi_1') ...
        ' -p ' fullfile(contourDir, 'prob_endo_epi_rv_rvi_1') ...
        ' --endo ' fullfile(contourDir, 'endo_epi_rv_rvi_endo_1') ... 
        ' --epi ' fullfile(contourDir, 'endo_epi_rv_rvi_epi_1') ... 
        ' --rv ' fullfile(contourDir, 'endo_epi_rv_rvi_rv_1') ... 
        ' -t endo_epi_rv_rvi -f ' endo_epi_rv_rvi_model ];
    command
    dos(command);

    command = ['gadgetron_CMR_ML_perf_seg_sax -i ' fullfile(trainingDataDir, 'Gd_resized_training_3') ... 
        ' -s ' fullfile(contourDir, 'seg_endo_epi_rv_rvi_2') ...
        ' -p ' fullfile(contourDir, 'prob_endo_epi_rv_rvi_2') ...
        ' --endo ' fullfile(contourDir, 'endo_epi_rv_rvi_endo_2') ... 
        ' --epi ' fullfile(contourDir, 'endo_epi_rv_rvi_epi_2') ... 
        ' --rv ' fullfile(contourDir, 'endo_epi_rv_rvi_rv_2') ... 
        ' -t endo_epi_rv_rvi -f ' endo_epi_rv_rvi_model ];
    command
    dos(command);
end
   
endo_epi_rv_rvi.rv0 = analyze75read(fullfile(contourDir, 'endo_epi_rv_rvi_rv_0'));
endo_epi_rv_rvi.rv1 = analyze75read(fullfile(contourDir, 'endo_epi_rv_rvi_rv_1'));
endo_epi_rv_rvi.rv2 = analyze75read(fullfile(contourDir, 'endo_epi_rv_rvi_rv_2'));

endo_epi_rv_rvi.prob_0 = analyze75read(fullfile(contourDir, 'prob_endo_epi_rv_rvi_0'));
endo_epi_rv_rvi.prob_1 = analyze75read(fullfile(contourDir, 'prob_endo_epi_rv_rvi_1'));
endo_epi_rv_rvi.prob_2 = analyze75read(fullfile(contourDir, 'prob_endo_epi_rv_rvi_2'));

endo_epi_rv_rvi.prob0_rv = endo_epi_rv_rvi.prob_0(:,:,4);
endo_epi_rv_rvi.prob1_rv = endo_epi_rv_rvi.prob_1(:,:,4);
endo_epi_rv_rvi.prob2_rv = endo_epi_rv_rvi.prob_2(:,:,4);

endo_epi_rv_rvi.prob0_rvi = endo_epi_rv_rvi.prob_0(:,:,5);
endo_epi_rv_rvi.prob1_rvi = endo_epi_rv_rvi.prob_1(:,:,5);
endo_epi_rv_rvi.prob2_rvi = endo_epi_rv_rvi.prob_2(:,:,5);

[maxP, I0] = max(endo_epi_rv_rvi.prob0_rvi(:));
[x0, y0] = ind2sub(size(endo_epi_rv_rvi.prob0_rvi), I0);
[maxP, I1] = max(endo_epi_rv_rvi.prob1_rvi(:));
[x1, y1] = ind2sub(size(endo_epi_rv_rvi.prob1_rvi), I1);
[maxP, I2] = max(endo_epi_rv_rvi.prob2_rvi(:));
[x2, y2] = ind2sub(size(endo_epi_rv_rvi.prob2_rvi), I2);

endo_epi_rv_rvi.rvi = [x0 y0; x1 y1; x2 y2];
