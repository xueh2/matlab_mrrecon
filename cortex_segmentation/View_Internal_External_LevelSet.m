
function [I, mapI] = View_Internal_External_LevelSet(I_internal, I_external, map, contour_intensity)
% combine the internal and external contours into one image

mapI = map;

% internal: red
% external: yellow

index2 = find(I_internal==contour_intensity);
I = I_external;

for i = 0:255
    index = find( (I_internal==i) | (I_external==i) );
    if( isempty(index) )
        break;
    end
end

I(index2) = i;
mapI(i+1, :) = [1 0 0];

imview(I, mapI);
return;