
function  D_pts = UseScalarField2D_BSpline2(pts_Image, Coeff)
% use the deformation
% SliceIndex starts from 0
% pts_Image starts from 0

num = size(pts_Image, 1);
numSlice = size(Coeff, 3);

D_pts = zeros(numSlice, num);

if ( ~isempty(Coeff) )
    for kk=1:numSlice
        for ii=1:num
            x = pts_Image(ii, 1);
            y = pts_Image(ii, 2);
            z = kk-1;

            D_pts(kk, ii) = Matlab_EvaluateBSplineInterpolation(Coeff, x, y, z, 'Row-wise');
        end
    end
end

return