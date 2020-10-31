
function [x, indexes] = kmean_init_MultiChannel(imagedata, header, brainmask, channelnumber)
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
nin = channelnumber;
x = zeros(ndata, nin); 
for ps = 1:ndata
    for m = 1:nin
        x(ps, m) = imagedata{m}(i(ps), j(ps), k(ps));
    end
end
indexes = [i j k];

return