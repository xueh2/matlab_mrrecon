
function [surface_points, num] = GetBoundaryPoints(data, level)
% get the boundary points (integer points)
% the 6 neighbors are used for determining the boundary points

[ysize, xsize, zsize] = size(data);

surface_points = zeros(floor(xsize*ysize*zsize/2), 3);
num = 0;

data = data - level;
index = find(data==0);
if ( isempty(index) == 0 )
    data(index) = -1e-4;
end

% header.xsize = xsize;
% header.ysize = ysize;
% header.zsize = zsize;
% 
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
% 
%  offsets6 =  [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
%  
% for k=2:zsize-1 % depth
%     for j=2:ysize-1 % row
%         for i=2:xsize-1 % col
%             
%             [NeighborStar, NeighborStar_Subs] = GetNeighborStar(j, i, k, [ysize xsize zsize], offsets6);
% %             neighborI = GetNeighbourhood(data, header, i, j, k, 3);
%             neighborI = data(NeighborStar);
%             index1 = find(neighborI>0);
%             index2 = find(neighborI<0);
%             
%             if ( isempty(index1)==0 )
%                 if ( isempty(index2)==0 )
%                     num = num+1;
%                     surface_points(num,:) = [i j k];
%                 end
%             end
%             
%         end
%     end
% end
% surface_points = surface_points(1:num,:);


% for x direction
for k=1:zsize
    for j=1:ysize
        for i=1:xsize-1
            y = j;
            z = k;
            % for points
            p1 = [i j k];
            v1 = data(j, i, k);
 
            p2 = [i+1 j k];
            v2 = data(j, i+1, k);

            if ( v1*v2<0 )
                if ( v1<0 )
                    num = num + 1;
                    surface_points(num,:) = p1;
                else
                    num = num + 1;
                    surface_points(num,:) = p2;
                end
            end
        end
    end
end

% for y direction
for k=1:zsize
    for i=1:xsize
        for j=1:ysize-1
            x = i;
            z = k;
            % for points
            p1 = [i j k];
            v1 = data(j, i, k);
 
            p2 = [i j+1 k];
            v2 = data(j+1, i, k);

            if ( v1*v2<0 )
                if ( v1<0 )
                    num = num + 1;
                    surface_points(num,:) = p1;
                else
                    num = num + 1;
                    surface_points(num,:) = p2;
                end
            end
        end
    end
end

% for z direction
for j=1:ysize
    for i=1:xsize
        for k=1:zsize-1
            x = i;
            y = j;
            % for points
            p1 = [i j k];
            v1 = data(j, i, k);
 
            p2 = [i j k+1];
            v2 = data(j, i, k+1);

            if ( v1*v2<0 )
                if ( v1<0 )
                    num = num + 1;
                    surface_points(num,:) = p1;
                else
                    num = num + 1;
                    surface_points(num,:) = p2;
                end
            end
        end
    end
end
surface_points = surface_points(1:num, :);
return