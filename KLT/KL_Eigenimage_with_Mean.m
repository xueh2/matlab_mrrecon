% This is to get the eigenimage only
% [a_eigen,V,D] = KL_Eigenimage_with_Mean(a);

function [a_eigen,V,D] = KL_Eigenimage_with_Mean(a);

%a = single(a);
s = size(a);
N = s(3);
if length(s) ~=3,
    'Input array is not a 3-D array!'
else
    a = (reshape(a, [s(1)*s(2),s(3)] )); %clear a
    A = (a'*a)/(s(1)*s(2)); %clear b;
    [V,D] = (eig(A)); % This is to diagonalize the covariance matrix
    V = (V);
    D = (D);
    a_eigen = (reshape( a*V , [s(1),s(2),s(3)] ) );
end










