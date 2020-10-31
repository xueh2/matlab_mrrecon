function [im, center, width] = autowindowing(im, xp, yp)
%--------------------------------------------------------------------------
% autowindowing: automatically detect the window size and center for best
% fit of image display.
%
%   [im, center, width] = autowindowing(im, xp, yp)
%
% Defined variables:
%   im (input):      image
%   xp:              x axis minimum percentage (window length)
%   yp:              y axis minimum percentage
%
%   im (output):     image modified to the range of 
%                    [center-width/2. center+width/2.]
%   center:          window center
%   width:           window width
%--------------------------------------------------------------------------

bins = 100;
step = (max(im(:))-min(im(:)))/bins;
t0 = numel(im);
%hist(im, bins);
[nelements, centers] = hist(im(:), bins);
s1 = struct('nelements',num2cell(nelements),'centers',num2cell(centers));
[~, order] = sort([s1(:).nelements],'descend');
s2 = s1(order);
t = 0;
for i = 1:numel(s2)
	t = t + s2(i).nelements;
    if (t / t0) > yp
        break;
	end
    if (max([s2(1:i).centers]) - min([s2(1:i).centers]))/(max(im(:))-min(im(:))) > xp
        break;
    end
end

if i == 1
    center = s2(1).centers;
    width = step;
else
    center = (max([s2(1:i).centers]) + min([s2(1:i).centers])) / 2.;
    width = (max([s2(1:i).centers]) - min([s2(1:i).centers])) + step;
end
im(im>(center+width/2.)) = (center+width/2.);
im(im<(center-width/2.)) = (center-width/2.);

end