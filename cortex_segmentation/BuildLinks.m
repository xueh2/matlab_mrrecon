
function links = BuildLinks(pts, numCells, cells)
% build the links as the vtkPolyData
% links is a struct array with the size of numOfPoints * 1
% everyElement = struct('ptID', 0, 'cellIDs', [])
% every element is a struct with two fields: ptID(scalar), cellIDs(N*1)
% cellIDs store the indexes to cells using the point
% the indexes mean the starting position for this cell in the 1D cells array

num = size(pts, 1);

links = cell(num, 1);

% oneElement = struct('ptID', 0, 'cellIDs', []);
IDs = pts(:,1);
for i=1:num
    links{i} = struct('ptID', 0, 'cellIDs', []);
    links{i}.ptID = IDs(i);
end

index = 1;

for i=1:numCells % for every cell
    cellLen = cells(index);
    
    cellIDs = cells(index+1:index+cellLen); % the pts used in current cell
    for k=1:cellLen
%         place = find(IDs==cellIDs(k));
%         if ( isempty(find(links{place}.cellIDs == i))==1 )
        place = cellIDs(k)+1;
        links{place}.cellIDs = [links{place}.cellIDs; index ];            
%         end
    end
    
    index = index + cellLen + 1;
end
return;