
function MaskVTK_PolyData_Binary(vtkfile, vtkfile_output, brainMask_filename)
% mask the vtkPolyData using the binary mask
% the criterion is the within rule

[mask, header] = LoadAnalyze(brainMask_filename, 'Grey');

[pts, numCells, cells] = VTKFile2PointCell(vtkfile);

pts_Image = pts;
pts_Image(:, 2:4) = world2image_DITK2(pts(:, 2:4), header);

pts_masked = MaskPts(mask, pts_Image);
pts_masked(:, 2:4) = image2world_DITK2(pts_masked(:, 2:4), header);

% ptIDs = pts_masked(:,1);
% save(matfilename{i}, 'ptIDs');

RefineCells_vtk(pts_masked, numCells, cells, vtkfile_output);
PolyConnectivity_vtk(vtkfile_output, vtkfile_output);
return;

