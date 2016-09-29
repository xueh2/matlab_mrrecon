function out=interpimages1(in,method,imsize);
% function out=interpimages1(in,method,imsize);
% function to interpolate image data
%
% Interpolation method:
%       'nearest' - nearest neighbor interpolation
%       'linear'  - bilinear interpolation (default)
%       'cubic'   - bicubic interpolation
%       'spline'  - spline interpolation
%
% Usage:
% 	out=interpimages1(in,method,imsize);
% 	out=interpimages1(in,method); % defaults to imsize=256


%      ***************************************
%      *  Peter Kellman  (kellman@nih.gov)   *
%      *  Laboratory for Cardiac Energetics  *
%      *  NIH NHLBI                          *
%      ***************************************

if nargin<2; method='linear'; imsize=256; end % set default method and imsize
if nargin<3; imsize=256; end % set default imsize

% check if input image size >= imsize
if size(in,1)>=imsize; return; end
if size(in,2)>=imsize; return; end

% determine aspect ratio (same for all axes)
aspect_ratio=size(in,1)/size(in,2);

% interpolate images
x=linspace(0,1,size(in,2));
y=linspace(0,1,size(in,1))';
if aspect_ratio >=1
    xi=linspace(0,1,floor(imsize/aspect_ratio));
    yi=linspace(0,1,imsize)';
else
    xi=linspace(0,1,imsize);
    yi=linspace(0,1,floor(imsize*aspect_ratio))';
end

for i=1:size(in,3)
    for j=1:size(in,4)
        for k=1:size(in,5)
            out(:,:,i,j,k)=interp2(x,y,in(:,:,i,j,k),xi,yi,method);
        end
    end
end

return