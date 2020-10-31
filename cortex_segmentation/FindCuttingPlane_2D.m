
function slicelabel = FindCuttingPlane_2D(slice, theta, transi)
% the left hemisphere is 1 and the right 2
% search range theta(1) ~ theta(2) (degree)
% transi(1) ~ transi(2), pixel
% the center is the image center

slice = double(slice);
slicelabel = zeros(size(slice), 'uint32');

[row, col] = size(slice);

lowTheta = theta(1);
highTheta = theta(2);

lowTransi = transi(1);
highTransi = transi(2);

corrR = zeros(length([lowTheta:highTheta]), 1);
transiR = zeros(length([lowTheta:highTheta]), 2); % (row, col), not(x, y)

max_corrR = -1;
max_transiR = [-1 -1];
index = 0;

reflected_slice = reflectSlice_2D(slice, 2); % 1: row reflection; 2: col reflection

for angle = lowTheta:highTheta
    angle
    index = index + 1;
    reflected_sliceRotated = imrotate(reflected_slice, angle, 'nearest', 'crop');
    
    for j = lowTransi:highTransi % col, horizontal direction
        for i = lowTransi:highTransi % row, vertical direction
            
            reflected_sliceRotated_Transi = imageTranslate(reflected_sliceRotated, j, i);
            
            % compute corrR
            current_corrR = ComputeDSC_2D(slice, reflected_sliceRotated_Transi);
            % select the maximal value
            if ( current_corrR > max_corrR )
                max_corrR = current_corrR;
                max_transiR = [i j];
            end
        end
    end
    
    corrR(index) = max_corrR;
    transiR(index, :) = max_transiR;
    max_corrR = -1;
    max_transiR = [-1 -1];
end

% compose the splitting line
[max_corrR, maxLabel] = max(corrR);
angle = lowTheta:highTheta;

optimal_angle = angle(maxLabel)/2; % degree
optimal_angle = optimal_angle*pi/180;
optimal_transi = transiR(maxLabel, :);

cosTheta = cos(optimal_angle);
sinTheta = sin(optimal_angle);
C = cosTheta*optimal_transi(2) + sinTheta*optimal_transi(1);

[rows, cols] = ind2sub(size(slice), find(slice>0));
num = length(rows);

for k = 1:num
    yI = rows(k);
    xI = cols(k);
    y = -(yI-1-(row-1)/2);
    x = xI-1-(col-1)/2;
    lc = cosTheta*x + sinTheta*y - C;
    if ( lc>=0 )
        slicelabel(yI,xI) = 1;
    else
        slicelabel(yI,xI) = 2;            
    end
end
slicelabel2 = slicelabel;
for tt = 1:2
    [row, col] = ind2sub(size(slice), find(slicelabel2==tt));
    meanCol(tt) = mean(col);
end

if ( meanCol(1)<meanCol(2) )
    slicelabel(find(slicelabel2==1)) = 1;
    slicelabel(find(slicelabel2==2)) = 2;
else
    slicelabel(find(slicelabel2==1)) = 2;
    slicelabel(find(slicelabel2==2)) = 1;
end

return;
