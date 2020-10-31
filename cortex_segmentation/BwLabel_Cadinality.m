
function T = BwLabel_Cadinality(cubicVolume, connectivity, adjacentFlag, Cubic_6, Cubic_18, Cubic_26)
% compute the topology number by counting the number of connected components (CC)
% adjacentFlag: 1 the effective CC must be adjacent to the center point; 
% 0 the effective CC can not be adjacent to the center point

ind = find(cubicVolume>0);
num = length(ind);
T = 0;

switch connectivity
    case 6
        neighbor = Cubic_6;
    case 18
        neighbor = Cubic_18;
    case 26
        neighbor = Cubic_26;
end

if ( num==0 )
    return;
end

if ( adjacentFlag==0 )
    [label, num] = bwlabeln(cubicVolume, connectivity);
    T = num;
end

if ( adjacentFlag==1 )
    [label, num] = bwlabeln(cubicVolume, connectivity);
    for tt=1:num
        index = find(label==tt);
        common_index = intersect(index, neighbor);
        if ( length(common_index)>0 )
            T = T+1;
        end
    end
end

return;