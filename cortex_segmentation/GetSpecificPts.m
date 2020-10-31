
function pts_Specific = GetSpecificPts(points, pt_IDs)
% points: 4*N
% pt_IDs : IDs, m ids
% pts_Specific: the 4*M points with the id in pt_IDs

m = length(pt_IDs);
pts_Specific = zeros(m, 4);
%pts_Specific(:,1) = pt_IDs;

pts_Specific = points(pt_IDs+1, :);

% points_IDs = points(:,1);
% place = 0;
% for i=1:m
%     index = find(points_IDs==pt_IDs(i));
%     if ( isempty(index) )
%         continue;
%     end
%     
%     place = place+1;
%     pts_Specific(place,:) = points(index, :);
% end
% pts_Specific = pts_Specific(1:place, :);
return;