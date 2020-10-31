
function pts_inflated = Inflation_OneIteration(pts, numCells, cells, links, lemda)
% perform the surface inflation using the relaxation operator
% links is a struct array with the size of numOfPoints * 1
% everyElement = struct('ptID', 0, 'cellIDs', [])

pts_inflated = pts;

num = size(pts, 1);
IDs = pts(:,1);

sumCellArea = 0;
centre = [0 0 0];

disp('inflation one iteration starts ...');

for i=1:num % for every point
    
    Ti = pts(i, 2:4);
    
    sumCellArea = 0;
    centre = [0 0 0];
    BbyC = [0 0 0];
    
    if ( pts(i,1) ~= links{i}.ptID )
        disp('some errors ...');
        continue;
    end
    
    cellIDs = links{i}.cellIDs;
    numCells = length(cellIDs);
    
    for k=1:numCells
        
        index = cellIDs(k);
        
%         cellLen = cells(index);
%         if ( cellLen == 3 ) % triangles
            
%             p1 = find(IDs==cells(index+1));
%             if ( isempty(p1)==1 )
%                 continue;
%             end
%             v1 = pts(p1,2:4);
% 
%             p2 = find(IDs==cells(index+2));
%             if ( isempty(p2)==1 )
%                 continue;
%             end
%             v2 = pts(p2,2:4);
% 
%             p3 = find(IDs==cells(index+3));
%             if ( isempty(p3)==1 )
%                 continue;
%             end
%             v3 = pts(p3,2:4);
            
            p1 = cells(index+1)+1;
            v1 = pts(p1,2:4);

            p2 = cells(index+2)+1;
            v2 = pts(p2,2:4);

            p3 = cells(index+3)+1;
            v3 = pts(p3,2:4);
            
            l1 = norm(v1-v2);
            l2 = norm(v2-v3);
            l3 = norm(v1-v3);

            Heron = (l1+l2+l3)/2;
            cellArea = sqrt(Heron*(Heron-l1)*(Heron-l2)*(Heron-l3));
            sumCellArea = sumCellArea + cellArea;
            
            centre = (v1+v2+v3)/3;
            BbyC = BbyC + cellArea.*centre;
%         end
    end
    
    averageTi = BbyC / (sumCellArea+eps);
    
    pts_inflated(i, 2:4) = (1-lemda).*Ti + lemda.*averageTi; 
end
disp('inflation one iteration ends ...');
disp(' ');
return;