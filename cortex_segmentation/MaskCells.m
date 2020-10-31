
function [pts_masked, numofCells_masked, cells_masked] = MaskCells(mask, pts, numofCells, cells)
% mask vtk cells

pts_masked = MaskPts(mask, pts);

[numofCells_masked, cells_masked] = RefineCells(pts_masked, numofCells, cells);

return;