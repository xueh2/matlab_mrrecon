
function [meanNorm, maxNorm] = analyzeDeformationField2DNoJacobian(dx, dy, header)

xsize = header.sizeX;
ysize = header.sizeY;
xvoxelsize = header.spacingX;
yvoxelsize = header.spacingY;

border = 1;

dx = dx*xvoxelsize;
dy = dy*yvoxelsize;

vecNorm = sqrt(dx.*dx+dy.*dy);
vecNorm2 = vecNorm(border+1:ysize-border, border+1:xsize-border);
meanNorm = mean(vecNorm2(:));
maxNorm = max(vecNorm2(:));
