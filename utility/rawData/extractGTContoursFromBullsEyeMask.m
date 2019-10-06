
function C = extractGTContoursFromBullsEyeMask(masks, num_sector)
% C = extractGTContoursFromBullsEyeMask(masks, num_sector)

for ii=1:num_sector
    t = zeros(size(masks));
    ind = find(masks==ii);
    t(ind) = 1;
    
    ratio = 2;
    
    t2 = imresize(t, ratio, 'bilinear');
    
    [starting_points, ending_points] = CCMS_Contour(t2, 0.65, 4, 0, [1 0 0], 2);
    
    starting_points = (starting_points-1)/ratio + 0.5;
    ending_points = (ending_points-1)/ratio + 0.5;
    
%     for tt=1:size(starting_points, 1)-1
%        line([starting_points(tt,1) ending_points(tt,1)],[starting_points(tt,2) ending_points(tt,2)], ...
%            'LineWidth', 1.0, 'Color', [1 0 0]); 
%     end
    
    C(ii) = struct('starting_points', starting_points, 'ending_points', ending_points);
end
