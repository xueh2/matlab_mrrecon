function x = softThresh(a,y,t)
% apply joint sparsity soft-thresholding

% number of coil
numOfCoil = size(y,4);
       
s = size(y);
% original LLL coefficient
LLL = y(1:s(1)/2, 1:s(2)/2, 1:s(3)/2, :);

absy = sqrt(sum(abs(y).^2, 4));
unity = y./(repmat(absy,[1,1,1, numOfCoil])+eps);

res = absy;

if ( a.temporalSoftThreshScaleRatio ~= 1 )
    % spatial direction
    res(:,:,1:s(3)/2, :) = absy(:,:,1:s(3)/2,:)-t;
    % temporal direction
    res(:,:,s(3)/2+1:s(3),:) = absy(:,:,s(3)/2+1:s(3),:)-a.temporalSoftThreshScaleRatio*t;
else
    res = absy-t;
end

res = (res + abs(res))/2;

x = unity.*repmat(res,[1,1,1,numOfCoil]);

x(1:s(1)/2, 1:s(2)/2, 1:s(3)/2, :) = LLL;
