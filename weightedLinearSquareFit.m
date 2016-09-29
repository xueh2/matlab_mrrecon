function [a, b] = weightedLinearSquareFit(index, I, weight)

% --------------------------------
% use the lscov

A = [index' ones(size(index')) ];
if ( size(I, 2) ~= 1 )
    D = I';
else
    D = I;
end

[x,se_x,mse] = lscov(A, D, abs(weight)'.*abs(weight)');

a = x(1);
b = x(2);