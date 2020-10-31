
function MappingPolyData_VTKFileMask(vtkfile, vtkfile_mask, Thres_MaskingDistance, filename)
% mask the vtkPolyData using the mask
% Thres_MaskingDistance is in mm;

[pts, numCells, cells] = VTKFile2PointCell(vtkfile);
[pts_mask, numCells_mask, cells_mask] = VTKFile2PointCell(vtkfile_mask);

[minDist, nearestPoints] = GetNearestPoints_VTK(pts(:,2:4), pts_mask(:,2:4));
index = find(minDist<=Thres_MaskingDistance);
pts_masked = pts(index(:), :);

ptIDs = pts_masked(:,1);
save MappingPolyData_VTKFileMask_ptIDs ptIDs

RefineCells_vtk(pts_masked, numCells, cells, filename);
PolyConnectivity_vtk(filename, filename);

return;
