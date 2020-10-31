
function data_noholes = statistical_filling(data, header, threshold)
% fix the small errors using the 3*3*3 statistical filter (Dale et al.,
% paper)

xsize = header.xsize;
ysize = header.ysize;
zsize = header.zsize;

Thres = floor((1-threshold)*27);
data_noholes = data;
for k = 2:zsize-1
    k
    for j = 2:ysize-1
        for i = 2:xsize-1
            
            label = data(j, i, k);
            
            neighbourhood = GetNeighbourhood(data, header, i, j, k, 3);
            
            numLabel = length(find(neighbourhood==label));
            
            if ( numLabel <= Thres )
                data_noholes(j, i, k) = 1 - label;
            end
        end
    end
end

return;