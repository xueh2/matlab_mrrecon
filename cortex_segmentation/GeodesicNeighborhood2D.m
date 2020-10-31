
function [GN_inds, GN_subs, cubicVolume] = GeodesicNeighborhood2D(i, j, B, objectLabel, connectivity, orderNum, ...
                                                                  offsets4, offsets8)
% compute the geodesic neighborhood of (i, j) with respect to B of order k
% connectivity: 4, 8

shape = size(B);
center_ind = sub2ind(shape, i, j);

[NeighborStar_8, NeighborStar_8_Subs] = GetNeighborStar2D(i, j, shape, offsets8);
NeighborStar_8_InterSect_X = InterSect_Topology_onlyInd(NeighborStar_8, B, objectLabel);
 
if ( connectivity == 4 )
    % k = 1;
    [NeighborStar_4, NeighborStar_4_Subs] = GetNeighborStar2D(i, j, shape, offsets4);
    GN_inds_k = InterSect_Topology_onlyInd(NeighborStar_4, B, objectLabel);
    
    GN_inds_temp = [];
    for tt=2:orderNum
        num = length(GN_inds_k);
        for ss = 1:num
            % for every y
            ind_y = GN_inds_k(ss);
            [y_i, y_j] = ind2sub(shape, ind_y);
            [y_NeighborStar_4, y_NeighborStar_4_Subs] = GetNeighborStar2D(y_i, y_j, shape, offsets4);
            y_Neighbor_4 = [y_NeighborStar_4; ind_y];
            y_InterSect = intersect(y_Neighbor_4, NeighborStar_8_InterSect_X);
            GN_inds_temp = union(GN_inds_temp, y_InterSect);
        end
        GN_inds_k = GN_inds_temp;
        GN_inds_temp = [];
    end
end

if ( connectivity == 8 )
    % k = 1;
    [NeighborStar_8, NeighborStar_8_Subs] = GetNeighborStar2D(i, j, shape, offsets8);
    GN_inds_k = InterSect_Topology_onlyInd(NeighborStar_8, B, objectLabel);
    
    GN_inds_temp = [];
    for tt=2:orderNum
        num = length(GN_inds_k);
        for ss = 1:num
            % for every y
            ind_y = GN_inds_k(ss);
            [y_i, y_j] = ind2sub(shape, ind_y);
            [y_NeighborStar_8, y_NeighborStar_8_Subs] = GetNeighborStar2D(y_i, y_j, shape, offsets8);
            y_Neighbor_8 = [y_NeighborStar_8; ind_y];
            y_InterSect = intersect(y_Neighbor_8, NeighborStar_8_InterSect_X);
            GN_inds_temp = union(GN_inds_temp, y_InterSect);
        end
        GN_inds_k = GN_inds_temp;
        GN_inds_temp = [];
    end
end

GN_inds = GN_inds_k;
num = length(GN_inds);
cubicVolume = zeros(3,3);
GN_subs = [];

if ( num>0 )
    [ri, rj] = ind2sub(shape, GN_inds);
    GN_subs = zeros(num, 2);
    GN_subs(:,1) = ri;
    GN_subs(:,2) = rj;
    
    Center_subs = ones(num, 2);
    Center_subs(:,1) = i;
    Center_subs(:,2) = j;

    currentSubs = GN_subs - Center_subs + 2;
    
    inds = sub2ind([3 3], currentSubs(:,1), currentSubs(:,2));
    cubicVolume(inds) = 1;
end
return;