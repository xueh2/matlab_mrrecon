
function [surface_points, num] = GetSurfacePoints(data, level)
% get the surface points using linear interpolation

[ysize, xsize, zsize] = size(data);

surface_points = zeros(floor(xsize*ysize*zsize/2), 3);
num = 0;

data = data - level;
index = find(data==0);
if ( isempty(index) == 0 )
    data(index) = -1e-4;
end

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
                x = LinearInterpolation(i, v1, i+1, v2, 0);
            
                num = num + 1;
                surface_points(num,:) = [x y z];
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
                y = LinearInterpolation(j, v1, j+1, v2, 0);
            
                num = num + 1;
                surface_points(num,:) = [x y z];
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
                z = LinearInterpolation(k, v1, k+1, v2, 0);
            
                num = num + 1;
                surface_points(num,:) = [x y z];
            end
        end
    end
end
surface_points = surface_points(1:num, :);
%-------------------------------------------------------------
function x = LinearInterpolation(x1, v1, x2, v2, v)
x = (x2*(v-v1)+x1*(v2-v))/(v2-v1);
%-------------------------------------------------------------