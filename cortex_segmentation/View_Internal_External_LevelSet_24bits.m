
function I = View_Internal_External_LevelSet_24bits(I_internal, I_external, contour_intensity)
% combine the internal and external contours into one image
% internal: red
% external: yellow

I = I_external;

width = size(I_internal, 2);
height = size(I_internal, 1);

for j = 1:height
    for i = 1:width
        colors = reshape(I_internal(j, i, :), [1 3]);
        p = (colors==contour_intensity);
        if ( sum(p) == 3 )
            I(j, i, :) = [255 0 0];
        end
        
        colors = reshape(I(j, i, :), [1 3]);
        p = (colors==contour_intensity);
        if ( sum(p) == 3 )
            I(j, i, :) = [0 0 255];
        end
    end
end

imview(I);
return;