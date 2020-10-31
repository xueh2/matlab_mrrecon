
function [ICI, GLN, ECI, MLN] = ComputeCortex_SecondOrder_Statistics(meanCurvature, gaussCurvature, header, ...
    pts_Image, numCells, cells, cellsArea)
% compute the second order statistics for cortical surface

meanCurvature = single(meanCurvature);
coefficient_H = col2row( ComputeCoefficients(meanCurvature) );

gaussCurvature = single(gaussCurvature);
coefficient_K = col2row( ComputeCoefficients(gaussCurvature) );

ICI = 0;
GLN = 0;
ECI = 0;
MLN = 0;

index = 1;
IDs = pts_Image(:,1);
for i=1:numCells
    cellLen = cells(index);
    
    if ( cellLen == 3 ) % trangles
        
        p1 = find(IDs==cells(index+1));
        if ( isempty(p1)==1 )
            continue;
        end
        v1 = pts_Image(p1,2:4);
        H_v1 = EvaluateSpline2(coefficient_H, 3, v1(1), v1(2), v1(3));
        K_v1 = EvaluateSpline2(coefficient_K, 3, v1(1), v1(2), v1(3));
        
        p2 = find(IDs==cells(index+2));
        if ( isempty(p2)==1 )
            continue;
        end
        v2 = pts_Image(p2,2:4);
        H_v2 = EvaluateSpline2(coefficient_H, 3, v2(1), v2(2), v2(3));
        K_v2 = EvaluateSpline2(coefficient_K, 3, v2(1), v2(2), v2(3));
        
        p3 = find(IDs==cells(index+3));
        if ( isempty(p3)==1 )
            continue;
        end
        v3 = pts_Image(p3,2:4);
        H_v3 = EvaluateSpline2(coefficient_H, 3, v3(1), v3(2), v3(3));
        K_v3 = EvaluateSpline2(coefficient_K, 3, v3(1), v3(2), v3(3));
        
        deltaArea = cellsArea(i);
        
        %--------------------------------------------------------------
        meanH = (H_v1+H_v2+H_v3)/3;
        meanK = (K_v1+K_v2+K_v3)/3;
        
        % ICI
        if ( meanK >= 0 )
            ICI = ICI + meanK*deltaArea;
        end
        
        % ECI
        ECI = ECI + 4*meanH*sqrt(abs(meanH*meanH-meanK))*deltaArea;
        
        % GLN
        GLN = GLN + meanK*meanK*deltaArea;
        
        % MLN
        MLN = MLN + meanH*meanH*deltaArea;
    end
    
    index = index + cellLen + 1;
end

GLN = sqrt(sum(cellsArea)*GLN);
MLN = sqrt(MLN);

return;
