
function [meanNorm, maxNorm, meanLogJac, maxLogJac, logJac] = analyzeDeformationField2D(dx, dy, header)

xsize = header.sizeX;
ysize = header.sizeY;
xvoxelsize = header.spacingX;
yvoxelsize = header.spacingY;

meanNorm = 0;
maxNorm = -1;
meanLogJac = 0;
maxLogJac = -1;

num = 0;

border = 1;

dx = dx*xvoxelsize;
dy = dy*yvoxelsize;

[dxdx, dxdy] = gradient(dx, xvoxelsize, yvoxelsize);
[dydx, dydy] = gradient(dy, xvoxelsize, yvoxelsize);

vecNorm = sqrt(dx.*dx+dy.*dy);
vecNorm2 = vecNorm(border+1:ysize-border, border+1:xsize-border);
meanNorm = mean(vecNorm2(:));
maxNorm = max(vecNorm2(:));

logJac = zeros(size(dx));

for j=border+1:ysize-border
    for i=border+1:xsize-border
        
        num = num + 1;
        
        vecX = dx(j, i);
        vecY = dy(j, i);
        
%         vecNorm = sqrt(vecX*vecX + vecY*vecY);
%         
%         meanNorm = meanNorm + vecNorm;
%         
%         if ( vecNorm > maxNorm )
%             maxNorm = vecNorm;
%         end
        
%         dxdx = (dx(j, i+1) - dx(j, i-1)) / (2*xvoxelsize);
%         dxdy = (dx(j+1, i) - dx(j-1, i)) / (2*yvoxelsize);
% 
%         dydx = (dy(j, i+1) - dy(j, i-1)) / (2*xvoxelsize);
%         dydy = (dy(j+1, i) - dy(j-1, i)) / (2*yvoxelsize);
        
        % vecJac = det([1.0+dxdx(j, i) dxdy(j, i); dydx(j, i) 1.0+dydy(j, i)]);
        vecJac = (1.0+dxdx(j, i))*(1.0+dydy(j, i)) - dxdy(j, i)*dydx(j, i);
        if( vecJac < eps) 
            vecJac = eps;
        end
        
        logVecJac = abs(log(vecJac));
        meanLogJac = meanLogJac + logVecJac;
        
        if ( logVecJac > maxLogJac )
            maxLogJac = logVecJac;
        end
        
        logJac(j, i) = logVecJac;
        
    end
end

% meanNorm = meanNorm/num;
meanLogJac = meanLogJac/num;
