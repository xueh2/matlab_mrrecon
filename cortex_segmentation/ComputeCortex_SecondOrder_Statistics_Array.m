
function [ICI_Array, GLN_Array, ECI_Array, MLN_Array] = ComputeCortex_SecondOrder_Statistics_Array(meanCurvature, gaussCurvature)
% compute the second order statistics for cortical surface

meanCurvature = single(meanCurvature);
coefficient_H = col2row( ComputeCoefficients(meanCurvature) );

gaussCurvature = single(gaussCurvature);
coefficient_K = col2row( ComputeCoefficients(gaussCurvature) );

ICI_Array = zeros(size(meanCurvature), 'single');
GLN_Array = zeros(size(meanCurvature), 'single');
ECI_Array = zeros(size(meanCurvature), 'single');
MLN_Array = zeros(size(meanCurvature), 'single');
        
% ICI
index = find(gaussCurvature>0);
ICI(index(:)) = gaussCurvature(index(:))/(4*pi);

% ECI
ECI = 4*meanCurvature*sqrt(abs(meanCurvature.*meanCurvature-gaussCurvature));

% GLN
GLN = abs(gaussCurvature);

% MLN
MLN = abs(meanCurvature);

return;
