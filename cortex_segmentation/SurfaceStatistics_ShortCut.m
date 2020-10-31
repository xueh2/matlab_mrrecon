
function [surface_area, CR, ICI, ECI, GLN, MLN] = SurfaceStatistics_ShortCut(vtkfile)

[surface_area, cellsArea] = GetSurfaceArea_vtk(vtkfile);
[pts, numOfCells, cells] = VTkFile2PointCell(vtkfile);
surface_area_CH = GetQHull(pts);
CR = surface_area / surface_area_CH;

[ICI, ECI, GLN, MLN, MeanC, GaussC] = Compute_SecondOrderStatistics_OutlierRemoval(vtkfile);