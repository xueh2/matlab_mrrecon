
function saveGTPlusExportImages(folderName, dataRole, data, pixelSize)
% save the gtplus create images

N = numel(size(data));

if ( nargin < 4 )
    pixelSize = ones(N, 1);
end

origin = zeros(N, 1);
header = CreateGtImageHeader(data, origin, pixelSize);

Matlab_gt_write_analyze(single(data), header, dataRole);