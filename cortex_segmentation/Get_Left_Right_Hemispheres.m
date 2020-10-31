
function [leftHemisphere, headerLeft, rightHemisphere, headerRight] = ...
    Get_Left_Right_Hemispheres(data, header)

leftHemisphere = zeros(size(data), 'uint32');
rightHemisphere = zeros(size(data), 'uint32');

leftHemisphere(find(data==1)) = 1;
rightHemisphere(find(data==2)) = 1;

% ================================================ %
label = [1];
xsize=header.xsize;
ysize=header.ysize;
zsize=header.zsize;

label3D = zeros(size(leftHemisphere), 'uint8');
num = length(label);
for i = 1:num
    label3D(find(leftHemisphere==label(i))) = 1;
end
largestComponent = zeros(size(label3D),'uint8');

[l,num] = bwlabeln(label3D,6);

volumes=zeros(num,1);
total = xsize*ysize*zsize;
for i=1:total
    if (l(i)>0)
        volumes(l(i))=volumes(l(i))+1;
    end
end
vd = sort(volumes, 'descend');
[largest,index] = max(volumes);

leftHemisphere(find(l~=index(1))) = 0;

[leftup, rightdown] = getBoundingBox_BinaryVolume(leftHemisphere);
[leftHemisphere, headerLeft]=getROI(leftHemisphere, header, leftup,rightdown);
% SaveAnalyze(ROI, headerROI, 'leftHemisphere.hdr', 'Grey');
% ================================================ %
label = [1];
xsize=header.xsize;
ysize=header.ysize;
zsize=header.zsize;

label3D = zeros(size(rightHemisphere), 'uint8');
num = length(label);
for i = 1:num
    label3D(find(rightHemisphere==label(i))) = 1;
end
largestComponent = zeros(size(label3D),'uint8');

[l,num] = bwlabeln(label3D,6);

volumes=zeros(num,1);
total = xsize*ysize*zsize;
for i=1:total
    if (l(i)>0)
        volumes(l(i))=volumes(l(i))+1;
    end
end
vd = sort(volumes, 'descend');
[largest,index] = max(volumes);

rightHemisphere(find(l~=index(1))) = 0;

[leftup, rightdown] = getBoundingBox_BinaryVolume(rightHemisphere);
[rightHemisphere, headerRight]=getROI(rightHemisphere, header, leftup,rightdown);
% SaveAnalyze(ROI, headerROI, 'rightHemisphere.hdr', 'Grey');

% ================================================ %
