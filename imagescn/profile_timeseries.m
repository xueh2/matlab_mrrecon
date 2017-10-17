function [out, profile_length]=profile_timeseries(I, profilematfilename, profile_no, image_no, pixel_spacing);
% function [out]=profile_timeseries(I, profilematfilename, profile_no, image_no, pixel_spacing);
%
% function to create a time series of image intensity profiles (along
% general x,y coordinates) defined by a profile matfile created by imagescn
% profile_length is in units of pixels, unless pixel_spacing is provided as
% input argument, and then it is in correct distance units.

%      ***************************************
%      *  Peter Kellman  (kellman@nih.gov)   *
%      *  Laboratory for Cardiac Energetics  *
%      *  NIH NHLBI                          *
%      ***************************************

if nargin <3
    profile_no=1; image_no=1;
end
if nargin < 5;
    pixel_spacing=1;
end

load (profilematfilename);

xpts=Profile_info_table(profile_no, image_no).Profile_Data_Px;
ypts=Profile_info_table(profile_no, image_no).Profile_Data_Py;
Profile_points=256;
if nargout==2
    profile_length=sqrt((xpts(1)-xpts(end))^2 + (ypts(1)-ypts(end))^2);
    profile_length=profile_length*pixel_spacing;
end

for i=1:size(I,3)
    [Px, Py, out(:,i)] = improfile(I(:,:,i),xpts,ypts,Profile_points,'bilinear');
end

return
