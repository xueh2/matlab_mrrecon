function  fa = dual_flip_angle_map(in, display_figure);
% function  fa = dual_flip_angle_map(in, display_figure);
%
% function to calculate the flip angle map (B1+ map)
% from dual flip angle images acquired at
% 2*theta and theta (e.g., 120 and 60 degrees)
%
% input:  in(row,col,FA_series)
%            where 1st series is FA = 2*theta, 2nd series is FA = theta
%         display_figure = 0 or 1, flag to display FA.
 
if nargin<2;
    display_figure = 1;
end
 
I=squeeze(in); % 
I=double(I);
 
for i=1:size(I,3)
    I(:,:,i)=medfilt2(I(:,:,i),[5 5]);
end
 
q=I(:,:,1)./(2*I(:,:,2));
q=clamp(q,[-1 1]);
fa=acos(q)*180/pi;
 
if display_figure
    figure; imagescn(fa); colormap('jet')
end
 
 
 
function out=clamp(in, value);
% function out=clamp(in, value);
%
% function will clamp the input to "value" for values > "value"
 
 
out=in;
 
if nargin<2
    return
elseif isempty(value)
    return
elseif length(value)==1 % assume value is maxvalue for clamp
    out(in > value)=value;
elseif length(value)==2 % assume value is vector of min and max values for clamp
    out(in > value(2)) = value(2);
    out(in < value(1)) = value(1);
end
