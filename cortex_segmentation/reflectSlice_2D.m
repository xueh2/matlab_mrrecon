
function reflected_slice = reflectSlice_2D(slice, dim)
% dim: 1, row reflection; 2, col reflection

reflected_slice = slice;

[row, col] = size(slice);

if ( dim == 2 )
    
    % col reflection
    for k = 1:col
        reflected_slice(:, k) = slice(:, col+1-k);
    end
end

if ( dim == 1 )
    
    % row reflection
    for k = 1:row
        reflected_slice(k, :) = slice(row+1-k, :);
    end
end