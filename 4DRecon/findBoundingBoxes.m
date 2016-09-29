

function [boundingBox, header] = findBoundingBoxes(volume, headers, pixelSpacing)
%
% this function finds the bounding box for a set of 2D spatial frames
%
% Inputs:
%    volume : a cell matrix, every element is a spatial volume (can be a 2D spatial frame)
%    headers : a cell array, header for every spatial volume
%    pixelSpacing : the pixel spacing for the header of bounding box
%
% Output:
%    boundingBox: [xmin xmax; ymin ymax; zmin zmax]
%    header: the header for boundingBox    
%
%     ***************************************
%     *  Hui Xue (hui-xue@siemens.com       *
%     *  2012-04                            *
%     ***************************************

N = numel(volume);

cp = zeros(8*N, 3);
for i=1:N
    [wc1, wc2, wc3, wc4, wc5, wc6, wc7, wc8] = computeEightCorner(volume{i}, headers{i});
    
    ind = (i-1)*8;
    cp(ind+1, :) = wc1;
    cp(ind+2, :) = wc2;
    cp(ind+3, :) = wc3;
    cp(ind+4, :) = wc4;
    cp(ind+5, :) = wc5;
    cp(ind+6, :) = wc6;
    cp(ind+7, :) = wc7;
    cp(ind+8, :) = wc8;    
end

xmin = min(cp(:,1));
xmax = max(cp(:,1));

ymin = min(cp(:,2));
ymax = max(cp(:,2));

zmin = min(cp(:,3));
zmax = max(cp(:,3));

boundingBox = zeros(3,2);
boundingBox(1,:) = [xmin xmax];
boundingBox(2,:) = [ymin ymax];
boundingBox(3,:) = [zmin zmax];

header = CreateFtkHeaderInfo(volume{1}, pixelSpacing);

header.positionPatient = [xmin ymin zmin];

rowVector = [xmax ymin zmin] - header.positionPatient;
rowVector = rowVector ./ norm(rowVector);

colVector = [xmin ymax zmin] - header.positionPatient;
colVector = colVector ./ norm(colVector);

header.orientationPatient = eye(3);
header.orientationPatient(1,:) = rowVector;
header.orientationPatient(2,:) = colVector;
header.orientationPatient(3,:) = cross(header.orientationPatient(1,:), header.orientationPatient(2,:));

header.sizeX = ceil( (xmax-xmin) / header.spacingX);
header.sizeY = ceil( (ymax-ymin) / header.spacingY);
header.sizeZ = ceil( (zmax-zmin) / header.spacingZ);
