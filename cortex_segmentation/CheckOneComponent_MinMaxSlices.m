
function labelHemsphere = CheckOneComponent_MinMaxSlices(labelHemsphere, headerLabel, minSlice, maxSlice, slicestep)
% check the connected component on every 2D slices.
% if one component have two different labels, the component is set to be
% the main label

zsize = headerLabel.zsize;

slices = [minSlice:minSlice+slicestep maxSlice-slicestep:maxSlice];
num = length(slices);
for k = 1:num
    slicelabel = slices(k);
    slice = labelHemsphere(:,:,slicelabel);
    
    sliceOld = double(slice>0);
    
    [l,num] = bwlabel(sliceOld, 8);

    if ( num <1 )
        continue;
    end
    
    for i = 1:num
        
        index = find(l==i);
        labels = slice(index(:));
        
        currentlabel = labels(1);
        
        index2 = find(labels~=currentlabel);
        if ( isempty(index2) == 0 )
            % have multiple labels
            area1 = length(find(labels==1));
            area2 = length(find(labels==2));
            
            if ( area1>area2 )
                slice(index(:)) = 1;
            else
                slice(index(:)) = 2;                
            end
        end
    end
    labelHemsphere(:,:,slicelabel) = slice;
end

return;