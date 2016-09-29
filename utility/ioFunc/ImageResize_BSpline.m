
function [dataZoomed, headerZoomed] = ImageResize_BSpline(data, scalingRatio)

S = size(data);
data2 = zeros(S(1), S(2), S(3)+4);
data2(:,:,1) = data(:,:,1);
data2(:,:,2) = data(:,:,1);
data2(:,:,3:end-2) = data;
data2(:,:,end-1) = data(:,:,end);
data2(:,:,end) = data(:,:,end);

header = CreateFtkHeaderInfo(data, [1 1 1 1]);
header2 = header;
header2.sizeZ = size(data2, 3);

headerVolume = header2;
headerVolume.sizeX = S(1) * scalingRatio;
headerVolume.sizeY = S(2) * scalingRatio;
headerVolume.sizeZ = S(3)+4;

headerVolume.spacingX = header.spacingX/scalingRatio;
headerVolume.spacingY = header.spacingY/scalingRatio;
dataVolume = zeros(headerVolume.sizeY, headerVolume.sizeX, headerVolume.sizeZ, 'single');

coeff = Matlab_ComputeBSplineCoefficient(double(data2), header2, 'Row-wise');        
[dataZoomed, headerZoomed] = Matlab_EvaluateBSplineResampledSlices(double(data2), coeff, header2, double(dataVolume), headerVolume, 'BSpline', 'Row-wise');
dataZoomed = dataZoomed(:, :, 3:end-2);
headerZoomed.sizeZ = size(dataZoomed, 3);
