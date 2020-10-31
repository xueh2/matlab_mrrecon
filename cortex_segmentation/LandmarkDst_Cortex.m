
function dists = LandmarkDst_Cortex(pts_World, surfaceVTKfile, levelsetFun, OriginalImageFile)

%dists = LandmarkDst_Cortex(pts_World, surfaceVTKfile, levelsetFun, OriginalImageFile);

[pts, numofcell, cells] = VTKFile2PointCell(surfaceVTKfile);

[dists, nearestPoints] = GetNearestPoints_VTK(pts_World, pts(:, 2:4));

return