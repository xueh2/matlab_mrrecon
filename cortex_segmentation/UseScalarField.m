
function  D_pts = UseScalarField(pts_Image, D, header)
% use the deformation

num = size(pts_Image, 1);

numSlice = header.zsize;
D_pts = zeros(num, numSlice);

for kk=1:numSlice
    % for every slice
    D_2D = D(:, :, kk);
    
    for tt=1:num
        % for every point
        x = pts_Image(tt, 1);
        y = pts_Image(tt, 2);
        
        local_d = interp2(D_2D, x+1, y+1, 'linear');
           
        D_pts(tt, kk) = local_d;
               
    end
end

return