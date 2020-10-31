
function [numCells_TH, cells_TH] = RefineCells(pts_TH, numCells, cells)
% remove a cell if it is using the points not belonging to pts_TH

numCells_TH = 0;
cells_TH = cells;
index_TH = 1;

ptIDs = pts_TH(:,1);

index = 1;
goodCell = true;
for k = 1:numCells
    
    cellLen = cells(index);
    
    goodCell = true;
    for j=1:cellLen
        if ( isempty(find(ptIDs==cells(index+j))) == 1 )
            goodCell = false;
            break;
        end
    end
    
    if ( goodCell == true )
        numCells_TH = numCells_TH + 1;
        cells_TH(index_TH) = cells(index);
        cells_TH(index_TH+1:index_TH+cellLen) = cells(index+1:index+cellLen);
        index_TH = index_TH + cellLen + 1;
    end
    
    index = index + cellLen + 1;
    
end

cells_TH = cells_TH(1:index_TH-1);

return;