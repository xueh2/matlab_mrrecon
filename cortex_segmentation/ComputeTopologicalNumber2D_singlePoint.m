
function [T_object, T_backgournd] = ComputeTopologicalNumber2D_singlePoint(i, j, B, connectivityObject, connectivityBackground)
% compute the topological numbers according to the algorithm proposed in 
% "Simple points, topological numbers and geodesic neighborhoods in cubic grids, PR Letters 15:1003-1011, 1994."
% B==1, object; B==0, background
% no boundary test, the levelset function must be within the effective data area
% for 2D digital topology, (connectivityObject, connectivityBackground) = {[4, 8], [8, 4]}

shape = size(B);

offsets4 =  [1 0 ; -1 0; 0 1; 0 -1];

offsets8 =  [-1    -1
     0    -1
     1    -1
    -1     0
     1     0
    -1     1
     0     1
     1     1];

if ( (connectivityObject==4) & (connectivityBackground==8) )
    % object, 4-connectivity
    [GN_inds, GN_subs, cubicVolume] = GeodesicNeighborhood2D(i, j, B, 1, connectivityObject, 2, offsets4, offsets8);
    T_object = BwLabel_Cadinality2D(cubicVolume, connectivityObject);

    % background, 8-connectivity
    [GN_inds, GN_subs, cubicVolume] = GeodesicNeighborhood2D(i, j, B, 0, connectivityBackground, 1, offsets4, offsets8);
    T_backgournd = BwLabel_Cadinality2D(cubicVolume, connectivityBackground);
    return;
end

if ( (connectivityObject==8) & (connectivityBackground==4) )      
    % object, 8-connectivity
    [GN_inds, GN_subs, cubicVolume] = GeodesicNeighborhood2D(i, j, B, 1, connectivityObject, 1, offsets4, offsets8);
    T_object(tt) = BwLabel_Cadinality2D(cubicVolume, connectivityObject);

    % background, 4-connectivity
    [GN_inds, GN_subs, cubicVolume] = GeodesicNeighborhood2D(i, j, B, 0, connectivityBackground, 2, offsets4, offsets8);
    T_backgournd(tt) = BwLabel_Cadinality2D(cubicVolume, connectivityBackground);
    return;
end

return;