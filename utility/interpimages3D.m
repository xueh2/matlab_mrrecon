function out=interpimages3D(in,method,imsize);
% function out=interpimages3D(in,method,[imsize1 imsize2 imsize3]);
% function to interpolate image data
%
% Interpolation method:
%       'nearest' - nearest neighbor interpolation
%       'linear'  - bilinear interpolation (default)
%       'cubic'   - bicubic interpolation
%       'spline'  - spline interpolation
%
% Usage:
% 	out=interpimages3D(in,method,imsize);

% interpolate images
x=linspace(0,1,size(in,2));
y=linspace(0,1,size(in,1))';
z=linspace(0,1,size(in,3))';

xi=linspace(0,1,imsize(2));
yi=linspace(0,1,imsize(1))';
zi=linspace(0,1,imsize(3))';


for j=1:size(in,4)
    for k=1:size(in,5)
        out(:,:,:,j,k)=interp3(x,y,z, in(:,:,:,j,k),xi,yi,zi,method);
    end
end

return