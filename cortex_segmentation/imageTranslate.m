
function reflected_sliceRotated_Transi = imageTranslate(reflected_sliceRotated, x_transi, y_transi)
% translate the image, centre coordinates
% no interpolation is needed

[row, col] = size(reflected_sliceRotated);
reflected_sliceRotated_Transi = zeros(size(reflected_sliceRotated));

for j = 1:col % target coordinates
    for i = 1:row
        
        tx = j - x_transi; % backmapping
        ty = i + y_transi;
        
        if ( (tx>0) & (tx<=col) & (ty>0) & (ty<=row) )
            reflected_sliceRotated_Transi(i, j) = reflected_sliceRotated(ty, tx);
        end
        
    end
end

