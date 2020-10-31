
function MaskingPolyData(vtkfile, mask, header, Thres_MaskingDistance, filename, matfilename)
% mask the vtkPolyData using the mask
% Thres_MaskingDistance is in mm;

[pts, numCells, cells] = VTKFile2PointCell(vtkfile);

pts_masked = MaskPts_DistanceThres(mask, header, pts, Thres_MaskingDistance);

ptIDs = pts_masked(:,1);
save(matfilename, 'ptIDs');

RefineCells_vtk(pts_masked, numCells, cells, filename);
PolyConnectivity_vtk(filename, filename);

return;
