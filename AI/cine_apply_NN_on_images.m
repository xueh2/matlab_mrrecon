
function cine_apply_NN_on_images(trainingDataDir, endo_model, epi_model, contourDir)

if(nargin<4)
    contourDir = trainingDataDir;
end
    
try
    tic
    Cine_resized_training = analyze75read(fullfile(trainingDataDir, 'Cine_resized_training_norm'));
    size(Cine_resized_training)
    toc
catch
    tic
    Cine_resized_training = analyze75read(fullfile(trainingDataDir, 'Cine_resized_training'));
    size(Gd_resized_training)
    toc
end

RO = size(Cine_resized_training, 1);
E1 = size(Cine_resized_training, 2);

PHS = size(Cine_resized_training, 3);
SLC = size(Cine_resized_training, 4);

% if(~isFileExist(fullfile(trainingDataDir, 'Cine_resized_training_1.img')))
    SLC = size(Cine_resized_training, 4);
    for slc=1:SLC
        header = CreateGtImageHeader(Cine_resized_training(:,:,:,slc));
        Matlab_gt_write_analyze(single(Cine_resized_training(:,:,:,slc)), header, fullfile(trainingDataDir, ['Cine_resized_training_' num2str(slc)]));
    end
% end

% return;
% ---------------------------------------------------------------------------------------------------------
% endo 

if(~isempty(endo_model))
    
    for slc=1:SLC
        command = ['gadgetron_CMR_ML_cine_seg_sax -i ' fullfile(trainingDataDir, ['Cine_resized_training_' num2str(slc)]) ... 
            ' -s ' fullfile(contourDir, ['seg_endo_' num2str(slc)]) ...
            ' -p ' fullfile(contourDir, ['prob_endo_' num2str(slc)]) ...
            ' --endo ' fullfile(contourDir, ['endo_' num2str(slc)]) ... 
            ' --endo-smoothing 48' ...
            ' -t endo -f ' endo_model ];
        command
        dos(command);
    end    
end

% ---------------------------------------------------------------------------------------------------------
% epi 
if(~isempty(epi_model))
    
    for slc=1:SLC
        command = ['gadgetron_CMR_ML_cine_seg_sax -i ' fullfile(trainingDataDir, ['Cine_resized_training_' num2str(slc)]) ... 
            ' -s ' fullfile(contourDir, ['seg_epi_' num2str(slc)]) ...
            ' -p ' fullfile(contourDir, ['prob_epi_' num2str(slc)]) ...
            ' --epi ' fullfile(contourDir, ['epi_' num2str(slc)]) ... 
            ' --endo-smoothing 48' ...
            ' -t epi -f ' epi_model ];
        command
        dos(command);
    end        
end
