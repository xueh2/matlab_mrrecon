
function [x, indexes] = kmean_init(imagedata, header, brainmask)
%kmean_init Initialises data vector and record the indexes
%
% Hui generated this as a part of neonatal cortex segmentation work

% new inputs: imagedata, header, brainmask, 
% the effective pixels are shown by 1 in the brainmask;
% find all effective points and record the [row col depth] into the indexes

% i, j, k: row, column, depth
[i, j, k] = ind2sub(size(brainmask), find(brainmask > 0));

points = [i, j, k];

i = points(:, 1);
j = points(:, 2);
k = points(:, 3);

ndata = length(i);
nin = 1; % single channel image
x = zeros(ndata, nin); 
for ps = 1:ndata
    x(ps, :) = imagedata(i(ps), j(ps), k(ps));
end
indexes = [i j k];

return