
function segResult = GetPVs_RectifyNonBrain(suspectWMPVs, segResult,...
    csflabel, cortexlabel, wmlabel, header);
% check the first-order neighbor of suspectWMPVs
% if the most neighboring voxels are gm, this volume is assigned as the incorrectly classified csf
% the prior of wm should decrease and the csf prior should increase

[l,num] = bwlabeln(suspectWMPVs,26);
sizeL = size(l);            

xsize = header.xsize;
ysize = header.ysize;
zsize = header.zsize;

disp(['Total region: ' num2str(num)]);

tic

seglabels = [];
indexInPrior = [];
offsets = [-1 0 0; 1 0 0; 0 -1 0; 0 1 0; 0 0 -1; 0 0 1];
numOffsets = size(offsets, 1);
disp(['Total region: ' num2str(num)]);
for ss = 1:num
    if ( mod(ss,50) == 0 )
        disp(['current region: ' num2str(ss)]);
    end
    [i, j, k] = ind2sub(sizeL, find(l==ss)); % row, col, depth
    numRegion = length(i);
    seglabels = [];
    for tt = 1:numRegion
        y = i(tt);
        x = j(tt);
        z = k(tt);
        index = 0;
        neighbors = zeros(numOffsets, 1);

        for pp = 1:numOffsets
            px = x + offsets(pp, 1); py = y + offsets(pp, 2); pz = z + offsets(pp, 3);                
            if ( pointInROI([px py pz], [1, 1, 1], [xsize ysize zsize]) == true )
                if (l(py, px, pz) > 0)
                    continue;
                end
                neighbors(index+1) = segResult(py, px, pz);
                index = index + 1;
            end
        end
        
        if ( index > 0 )
            seglabels = [seglabels; neighbors(1:index)];
        end
    end
    if ( isempty(seglabels) == 0 )
       
        recordflag = false;
        
        numGM = length(find(seglabels==cortexlabel));
        numCSF = length(find(seglabels==csflabel));
        numWM = 0;
        lenWM = length(wmlabel);
        for tt = 1:lenWM
            numWM = numWM + length(find(seglabels==wmlabel(tt)));
        end

        [maxTissue, I] = max([numGM, numCSF, numWM]);
        if ( I == 1 )
            maxL = cortexlabel;
        end
        
        if ( I == 2 )
            maxL = csflabel;
        end

        if ( I == 3 )
            maxL = wmlabel(1);
        end
        
        segResult(find(l==ss)) = maxL;
        
    end
end
toc
return;