function [y] = function_Phi2D(x1,x2,m)

% function of polynom from 2 variables phi(x1,x2) of order m
%
%   [y] = function_Phi2D(x1,x2,m)
%
%   x1, x2 - are the scalar arguments of function phi
%
%   m - is an order of polynom
% 
% Dmitriy Paliy. Tampere University of Technology. TICSP. 16-02-2005
% dmitriy.paliy@tut.fi

switch m
    case 1
        y = 1;
    case 2
        y = x1;
    case 3
        y = x2;
    case 4
        y = x1^2/2;
    case 5
        y = x2^2/2;
    case 6
        y = x1*x2;
    case 7
        y = x1^3/6;
    case 8
        y = x2^3/6;
    case 9
        y = x2*x1^2/2;
    case 10
        y = x1*x2^2/2;
end;    