function [points, F] = plotFlowMask(mag, mask, points)

PHS = size(mag, 3);

if nargin < 3
    points = cell(PHS, 2);

    for ii=1:PHS
        [starting_points, ending_points] = CCMS_Contour(mask(:,:,ii), 0.5, 8, 0, [1 0 0], 1);    
        points{ii,1} = starting_points;
        points{ii,2} = ending_points;
    end
end

f = figure;
for ii=1:PHS
    f = RenderContour(mag(:,:,ii), points{ii,1}, points{ii,2}, [1 0 0], 1, f);
    F(ii) = getframe(f);
end
