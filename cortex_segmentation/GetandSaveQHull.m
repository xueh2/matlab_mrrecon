
function surface_area_CH = GetandSaveQHull(vtkfile, QHullFileName)
% compute the convex hull
% pts is in the world coordinate (mm)

[pts, numCells, cells] = VTKFile2PointCell(vtkfile);

X = pts(:,2:4);
C = convhulln(X);

num = size(C,1);
surface_area_CH = 0;

for i=1:num
    
    p1 = C(i,1);
    v1 = pts(p1,2:4);
        
    p2 = C(i,2);
    v2 = pts(p2,2:4);
    
    p3 = C(i,3);
    v3 = pts(p3,2:4);

    l1 = norm(v1-v2);
    l2 = norm(v2-v3);
    l3 = norm(v1-v3);

    Heron = (l1+l2+l3)/2;
    surface_area_CH = surface_area_CH + sqrt(Heron*(Heron-l1)*(Heron-l2)*(Heron-l3));
end

% save QHull as a vtk triangle mesh
newPts = zeros(num*3, 4);
newPts(:, 1) = 0:(num*3-1);

numCells = num;
newCells = zeros(4*num, 1, 'int32');

index = 1;
for i=1:num
   
    p1 = C(i,1);
    v1 = pts(p1,2:4);
        
    p2 = C(i,2);
    v2 = pts(p2,2:4);
    
    p3 = C(i,3);
    v3 = pts(p3,2:4); 
    
    newPts(3*(i-1)+1, 2:4) = v1;
    newPts(3*(i-1)+2, 2:4) = v2;
    newPts(3*(i-1)+3, 2:4) = v3;
    
    newCells(4*(i-1)+1) = 3;
    newCells(4*(i-1)+2) = 3*(i-1);
    newCells(4*(i-1)+3) = 3*(i-1)+1;
    newCells(4*(i-1)+4) = 3*(i-1)+2;
end

PointCell2VTKFile(newPts, num, int32(newCells), QHullFileName);
PolyDataClean(QHullFileName, QHullFileName);
return;