
function [ICI, GLN, ECI, MLN] = ComputeCortex_SecondOrder_Statistics_VTK(meanC, gaussC, ...
    pts, numCells, cells, cellsArea)
% compute the second order statistics for cortical surface
% the meanC and gaussC have been computed using the vtkCurvatures

ICI = 0;
GLN = 0;
ECI = 0;
MLN = 0;

index = 1;
IDs = pts(:,1);
for i=1:numCells
    cellLen = cells(index);
    
    if ( cellLen == 3 ) % trangles
        
        p1 = find(IDs==cells(index+1));
        if ( isempty(p1)==1 )
            continue;
        end
        v1 = pts(p1,2:4);
        H_v1 = meanC(p1);
        K_v1 = gaussC(p1);
        
        p2 = find(IDs==cells(index+2));
        if ( isempty(p2)==1 )
            continue;
        end
        v2 = pts(p2,2:4);
        H_v2 = meanC(p2);
        K_v2 = gaussC(p2);
        
        p3 = find(IDs==cells(index+3));
        if ( isempty(p3)==1 )
            continue;
        end
        v3 = pts(p3,2:4);
        H_v3 = meanC(p3);
        K_v3 = gaussC(p3);
        
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
