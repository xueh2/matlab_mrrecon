
function [ yOut, schemeDataOut ] = GetIntermediateDigitalTopology(t, yIn, schemeDataIn)
% schemeDataIn.LastData
% schemeDataIn.LastY
% schemeDataIn.B
% schemeData.connectivityObject
% schemeData.connectivityBackground

disp(['current time is ' num2str(t)]);
disp(['checking topology 2D ... ']);

% yOut phi at next time point
% yIn phi_temp
yOut = yIn;
schemeDataOut = schemeDataIn;

B = schemeDataIn.B; % B will be improved at every iteration
LastData = schemeDataIn.LastData;
tempData = reshape(yIn, schemeDataIn.shape);
connectivityObject = schemeDataIn.connectivityObject;
connectivityBackground = schemeDataIn.connectivityBackground;

global ntimes

map = gray(256);
simplePoints = zeros(size(B), 'uint8');
simplePoints(find(B==1)) = 80;
filename = ['./simplePoints2D/B_' num2str(ntimes) '.bmp'];
imwrite(simplePoints, map, filename, 'bmp');

yOut = yIn;

NextData = reshape(yOut, schemeDataOut.shape);

schemeDataOut.LastData = NextData;
schemeDataOut.LastY = yOut;
schemeDataOut.B = B;

schemeDataOut = processTopology(schemeDataOut);

ntimes = ntimes+1;
return;
