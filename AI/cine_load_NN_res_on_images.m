
function [endo, epi] = cine_load_NN_res_on_images(trainingDataDir, contourDir, load_endo, load_epi)

if(nargin<2)
    contourDir = trainingDataDir;
end

if(nargin<3)
    load_endo = 1;
end

if(nargin<4)
    load_epi = 1;
end
    
endo = struct('seg', [], 'prob', [], 'C', []);
epi = struct('seg', [], 'prob', [], 'C', []);

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

endo.seg = zeros(RO, E1, PHS, SLC);
endo.prob = zeros(RO, E1, PHS, SLC);

for slc=1:SLC
    endo.seg(:,:,:,slc) = analyze75read(fullfile(contourDir, ['seg_endo_' num2str(slc)]));
    endo.prob(:,:,:,slc) = analyze75read(fullfile(contourDir, ['prob_endo_' num2str(slc)]));
end

% ---------------------------------------------------------------------------------------------------------
% endo 

if(load_endo)
    endo.C = cell(SLC, PHS);
    
    for slc=1:SLC
        for phs=1:PHS
            seg = endo.seg(:,:,phs,slc);
            if(max(seg(:))>0.5)
                try
                    C = mask2contour(seg, 1, 400, 48);    
                    endo.C{slc, phs} = C;
                catch
                    endo.C{slc, phs} = [];
                end
            end
        end
    end
end

% ---------------------------------------------------------------------------------------------------------
% epi 
if(load_epi)
       
    epi.seg = zeros(RO, E1, PHS, SLC);
    epi.prob = zeros(RO, E1, PHS, SLC);
    
    for slc=1:SLC
        epi.seg(:,:,:,slc) = analyze75read(fullfile(contourDir, ['seg_epi_' num2str(slc)]));
        epi.prob(:,:,:,slc) = analyze75read(fullfile(contourDir, ['prob_epi_' num2str(slc)]));
        
        if(load_endo)
            endo_seg = endo.seg(:,:,:,slc);
            epi_seg = epi.seg(:,:,:,slc);
            
            ind = find(endo_seg(:)>0);
            ind2 = find(epi_seg(:)>0);
            
            if(numel(ind2)>numel(ind))
                epi_seg(ind) = 1;
                epi.seg(:,:,:,slc) = epi_seg;
            end
            
            if(numel(ind)>numel(ind2))
                ind3 = find(endo_seg(:)>0 & epi_seg(:)==0);
                endo_seg(ind3) = 0;
                endo.seg(:,:,:,slc) = endo_seg;
            end
        end
    end
    
    epi.C = cell(SLC, PHS);
    
    for slc=1:SLC
        for phs=1:PHS
            seg = epi.seg(:,:,phs,slc);
            if(max(seg(:))>0.5)
                try
                    C = mask2contour(seg, 1, 400, 48);    
                    epi.C{slc, phs} = C;
                catch
                    epi.C{slc, phs} = [];
                end
            end
        end
    end
end
