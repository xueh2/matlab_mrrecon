
function dists = LandmarkDst_Cortex_SDF(pts_World, SDF_File)
% compute the distance using the SDF function
% the negative distance means the point is inside the surface; positive means outside

p = dir(SDF_File);
if ( isempty(p)==1 )
    error('Can not find SDF_File ...');
    return;
end

[data, header] = LoadAnalyze(SDF_File, 'Real');

pts_Image = world2image_DITK2(pts_World, header);

SDF = single(data);
coefficient = col2row( ComputeCoefficients(SDF) );

num = size(pts_Image, 1);
dists = zeros(num, 1);

for i = 1:num
    x = pts_Image(i, 1);
    y = pts_Image(i, 2);
    z = pts_Image(i, 3);
    dists(i) = EvaluateSpline2(coefficient, 3, x, y, z);
end

voxelsize = (header.xvoxelsize + header.yvoxelsize + header.zvoxelsize)/3;
dists = dists .* voxelsize;

return;
