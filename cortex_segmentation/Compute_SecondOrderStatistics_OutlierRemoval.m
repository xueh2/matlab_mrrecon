
function [ICI, ECI, GLN, MLN, newMeanC, newGaussC] = Compute_SecondOrderStatistics_OutlierRemoval(vtkfile)

% [ptIDs, meanC, gaussC] = Curvature_PolyData(vtkfile);
% 
% N = 4096;
% minBinNumber = length(meanC)/N;
% if ( minBinNumber>=10 )
%     minBinNumber = 10;
% end
% 
% [low_meanC, high_meanC] = findThreshold_Histogram(meanC, minBinNumber, N)
% [low_gaussC, high_gaussC] = findThreshold_Histogram(gaussC, minBinNumber, N)
% 
% [ICI, ECI, GLN, MLN, newMeanC, newGaussC] = SecondOrder_Statistics_vtk_Threshold(vtkfile, low_meanC, high_meanC, low_gaussC, high_gaussC);

[ICI, ECI, GLN, MLN, newMeanC, newGaussC] = SecondOrder_Statistics_vtk(vtkfile);


return;