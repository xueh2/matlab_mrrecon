
function [surface_area, cellsArea] = GetSurfaceArea_vtkPloyData(pts, numCells, cells)
% compute the surface area of the vtkPloyData surfaces
% the pts is in the world coordinate (mm)

surface_area = 0;
cellsArea = zeros(numCells, 1);

index = 1;
IDs = pts(:,1);
for i=1:numCells
    cellLen = cells(index);
    
    if ( cellLen == 3 ) % trangles
        
        p1 = find(IDs==cells(index+1));
        if ( isempty(p1)==1 )
            continue;
        end
        v1 = pts(p1,2:4);
        
        p2 = find(IDs==cells(index+2));
        if ( isempty(p2)==1 )
            continue;
        end
        v2 = pts(p2,2:4);
        
        p3 = find(IDs==cells(index+3));
        if ( isempty(p3)==1 )
            continue;
        end
        v3 = pts(p3,2:4);
        
        l1 = norm(v1-v2);
        l2 = norm(v2-v3);
        l3 = norm(v1-v3);
        
        Heron = (l1+l2+l3)/2;
        cellsArea(i) = sqrt(Heron*(Heron-l1)*(Heron-l2)*(Heron-l3));
        surface_area = surface_area + cellsArea(i);
    end
    
    index = index + cellLen + 1;
end

return;