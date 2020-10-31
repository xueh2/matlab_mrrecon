
function [T_object, T_backgournd] = ComputeTopologicalNumber_singlePoint(i, j, k, B, ...
    connectivityObject, connectivityBackground, ...
    offsets6, Cubic_6, offsets18, Cubic_18, offsets26, Cubic_26)
% compute the topological numbers according to the algorithm proposed in 
% "Simple points, topological numbers and geodesic neighborhoods in cubic grids, PR Letters 15:1003-1011, 1994."
% B==1, object; B==0, background
% no boundary test, the levelset function must be within the effective data area
% for 3D digital topology, (connectivityObject, connectivityBackground) = {[6, 26], [26, 6], [6+, 18], [18, 6+]}

shape = size(B);

% offsets6 =  [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
% [NeighborStar_6, NeighborStar_6_Subs] = GetNeighborStar(2, 2, 2, [3 3 3], offsets6);
% Cubic_6 = [NeighborStar_6; 14]; % center point
% 
% offsets18 =  [0    -1    -1
%      -1     0    -1
%      0     0    -1
%      1     0    -1
%      0     1    -1
%     -1    -1     0
%      0    -1     0
%      1    -1     0
%     -1     0     0
%      1     0     0
%     -1     1     0
%      0     1     0
%      1     1     0
%      0    -1     1
%     -1     0     1
%      0     0     1
%      1     0     1
%      0     1     1];
% [NeighborStar_18, NeighborStar_18_Subs] = GetNeighborStar(2, 2, 2, [3 3 3], offsets18);
% Cubic_18 = [NeighborStar_18; 14]; % center point
% 
% offsets26 =  [-1    -1    -1
%      0    -1    -1
%      1    -1    -1
%     -1     0    -1
%      0     0    -1
%      1     0    -1
%     -1     1    -1
%      0     1    -1
%      1     1    -1
%     -1    -1     0
%      0    -1     0
%      1    -1     0
%     -1     0     0
%      1     0     0
%     -1     1     0
%      0     1     0
%      1     1     0
%     -1    -1     1
%      0    -1     1
%      1    -1     1
%     -1     0     1
%      0     0     1
%      1     0     1
%     -1     1     1
%      0     1     1
%      1     1     1];
% [NeighborStar_26, NeighborStar_26_Subs] = GetNeighborStar(2, 2, 2, [3 3 3], offsets26);
% Cubic_26 = [NeighborStar_26; 14]; % center point

% [i_subs, j_subs, k_subs] = ind2sub(shape, index);
% num = length(index);
% T_object = ones(num, 1);
% T_backgournd = ones(num, 1);

if ( (connectivityObject==26) & (connectivityBackground==6) )
        
    [NeighborStar_26, NeighborStar_26_Subs] = GetNeighborStar(i, j, k, shape, offsets26);
    [InterResult, cubicVolume] = InterSect_Topology(NeighborStar_26, NeighborStar_26_Subs, i, j, k, B, 1); % object
    T_object = BwLabel_Cadinality(cubicVolume, connectivityObject, 0, Cubic_6, Cubic_18, Cubic_26);

    [NeighborStar_18, NeighborStar_18_Subs] = GetNeighborStar(i, j, k, shape, offsets18);
    [InterResult, cubicVolume] = InterSect_Topology(NeighborStar_18, NeighborStar_18_Subs, i, j, k, B, 0); % object
    T_backgournd = BwLabel_Cadinality(cubicVolume, connectivityBackground, 1, Cubic_6, Cubic_18, Cubic_26);
    return;
end

if ( (connectivityObject==6) & (connectivityBackground==26) )
    disp('not implemented yet ...');
    return;
end

if ( (connectivityObject==18) & (connectivityBackground==6) )
    % object, 18-connectivity
    [GN_inds, GN_subs, cubicVolume] = GeodesicNeighborhood3D(i, j, k, B, 1, connectivityObject, 2, offsets6, offsets18, offsets26);
%     T_object = BwLabel_Cadinality(cubicVolume, connectivityObject, 0, Cubic_6, Cubic_18, Cubic_26);
    T_object = BwLabel_Cadinality3D(cubicVolume, connectivityObject);

    % background, 6+ connectivity
    [GN_inds, GN_subs, cubicVolume] = GeodesicNeighborhood3D(i, j, k, B, 0, connectivityBackground, 3, offsets6, offsets18, offsets26);
%     T_backgournd = BwLabel_Cadinality(cubicVolume, connectivityBackground, 0, Cubic_6, Cubic_18, Cubic_26);
    T_backgournd = BwLabel_Cadinality3D(cubicVolume, connectivityBackground);
    return;
end

if ( (connectivityObject==6) & (connectivityBackground==18) )
    disp('not implemented yet ...');
    return;
end

return;