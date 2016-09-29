
function header = CreateGtImageHeader(data, origin, pixelSize, axis)
% header = CreateGtImageHeader(data, origin, pixelSize, axis)
% header = CreateGtImageHeader(data)

s = size(data);
nDim = numel(s);

if ( nargin < 2 )
    origin = zeros([1 nDim]);
end

if ( nargin < 3 )
    pixelSize = ones([1 nDim]);
end

if ( nargin < 4 )
    axis = eye(nDim, nDim);
end

header = struct('origin', single(origin), 'pixelSize', single(pixelSize), 'axis', single(axis));
