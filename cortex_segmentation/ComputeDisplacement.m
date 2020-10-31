
function  [Dx, Dy, Dz] = ComputeDisplacement(target_FullName, dofname, InverseFlag)
% compute the displacement of every voxel under a DITK transformation
% the displacement is in the voxel coordinate
% NewPoint_inImage = OldPoint_inImage + [Dx Dy Dz]

Dx = [];
Dy = [];
Dz = [];

if ( isempty(dir(target_FullName)) )
    disp(['Not existing : ' target_FullName]);
    return
end

if ( isempty(dir(dofname)) )
    disp(['Not existing : ' dofname]);
    return
end

[data, header] = LoadAnalyze(target_FullName, 'Grey');

Dx = zeros(size(data));
Dy = zeros(size(data));
Dz = zeros(size(data));

pts = zeros(header.xsize*header.ysize*header.zsize, 3);


index = 1;
for k=1:header.zsize
    for j=1:header.ysize
        for i=1:header.xsize
           
            pts(index, :) = [i-1 j-1 k-1];
            index = index + 1;
        end
    end
end

pts = Image2World_DITK2(pts, header);
pts(:, 3) = 0;
pts_Transformed = pointsTransformation(dofname, pts, InverseFlag);
pts_Transformed = World2Image_DITK2(pts_Transformed, header);

index = 1;
for k=1:header.zsize
    for j=1:header.ysize
        for i=1:header.xsize
            
            Dx(j, i, k) = pts_Transformed(index, 1) - i + 1;
            Dy(j, i, k) = pts_Transformed(index, 2) - j + 1;
            Dz(j, i, k) = pts_Transformed(index, 3) - k+ 1;
            
            index = index + 1; 
        end
    end
end

return