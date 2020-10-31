
function  D_pts = UseScalarField2D_BSpline(pts_Image, Coeff_Dx, Coeff_Dy, SliceIndex)
% use the deformation
% SliceIndex starts from 0
% pts_Image starts from 0

num = size(pts_Image, 1);
D_pts = pts_Image;

if ( ~isempty(Coeff_Dx) )
    for kk=1:num
        x = pts_Image(kk, 1);
        y = pts_Image(kk, 2);
        z = SliceIndex;
        
        D_pts(kk, 1) = x + EvaluateSpline2(Coeff_Dx, 3, x, y, z);
    end
end
if ( ~isempty(Coeff_Dy) )
    for kk=1:num
        x = pts_Image(kk, 1);
        y = pts_Image(kk, 2);
        z = SliceIndex;
        
        D_pts(kk, 2) = y + EvaluateSpline2(Coeff_Dy, 3, x, y, z);
    end
end

return