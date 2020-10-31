

function [Label_Hemsphere, minSlice, maxSlice] = LabelHemisphere(data, header, theta, transi)
% the left hemisphere is label as 1 and the right as 2;

zsize = header.zsize;

Label_Hemsphere = zeros(size(data), 'uint8');
meanR = [0 0];

minSlice = header.zsize;
maxSlice = -1;
for k = 1:zsize
% for k = 79:83
    
    slice = data(:,:,k);
    slicelabel = zeros(size(slice), 'uint8');
    
    [l,num] = bwlabel(slice, 8);

    if ( num > 2 )
        slice2 = bwareaopen(slice, 15);
        [l,num] = bwlabel(slice2, 8);
    end
    
    if ( num<1 )
        continue;
    end
    
    if ( k<minSlice )
        minSlice = k;
    end
    if ( k>maxSlice )
        maxSlice = k;
    end

    if ( num == 2 )
        for tt = 1:num
            [row, col] = ind2sub(size(slice), find(l==tt));
            meanCol(tt) = mean(col);
        end
        
        if ( meanCol(1)<meanCol(2) )
            slicelabel(find(l==1)) = 1;
            slicelabel(find(l==2)) = 2;
        else
            slicelabel(find(l==1)) = 2;
            slicelabel(find(l==2)) = 1;
        end
        Label_Hemsphere(:,:,k) = slicelabel;
    else
        disp(['current slice is ' num2str(k) ' ... ']);
        % perform a complex search
        slicelabel = FindCuttingPlane_2D(slice, theta, transi);
        Label_Hemsphere(:,:,k) = slicelabel;
    end
end
return;