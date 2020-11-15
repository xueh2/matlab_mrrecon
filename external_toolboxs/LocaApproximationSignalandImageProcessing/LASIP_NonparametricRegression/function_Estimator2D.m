function [yh] = function_Estimator2D(z,x,h,m,WindowType,derivative);

% The LPA estimate of m-th order for the 2D signal z(X_s) on the 
% IRREGULAR grid X_s.
%
% SYNTAX
%   yh = function_Estimator2D(z,x,h,m,WindowType,derivative);
%   yh = function_Estimator2D(z,x,h,m,WindowType);
%   yh = function_Estimator2D(z,x,h,m);
%   yh = function_Estimator2D(z,x,h);
% 
% DESCRIPTION
%   The function returns the LPA estimate of a signal or its derivative. 
%   The grid of arguments X_s where observations z are made is irregular in 
%   general. The grid of arguments x where the estimate yh s done is 
%   irregular as well. The estimate is obtained by the order m 
%   approximation with a scale h.
%
%   Here:
%     z - is a 3xn array consists of:
%       a) z(1,:) = X1_s; is a grid of arguments X1_s where z is observed. 
%          It has length n and can be regular or irregular.
%       b) z(2,:) = X2_s; is a grid of arguments X2_s where z is observed. 
%          It has length n and can be regular or irregular.
%       c) z(3,:) = z; function values.
% 
%     x - is a CELL array corresponds to an array of ARGUMENTS where the 
%         LPA estimation is done and, in general, is different from X_s 
%         and also irregular.
%     x{1} = x1 and x{2} = x2 arrays which define rectangular grid, but
%         they can have different length.
% 
%     h - estimator bandwidth (any real);
%
%     m - is an order of estimator bandwidth (m=0,1,2,3), m = 0 by default.
%
%     WindowType: 'Gaussian' (default), 'GaussianLeft', 'GaussianRight',
%                 'Rectangular', 'RectangularLeft', 'RectangularRight',
%                 'Treecube', 'Hermitian',
%                 'Exponential', 'ExponentialLeft', 'ExponentialRight',
%                 'InterpolationWindow1', 'InterpolationWindow2', 
%                 'InterpolationWindow3'
%
%     derivative - (r in the book notations) is an order of estimation 
%         derivative, which is 0 by default.
% 
% RETURNS
%   yh - is a 2D estimated signal
% 
% REMARKS
%   The example is illustarated on Figure 2.11 of the book. For more details
%   read 2.2.5 Examples part p. 41.
%
% Dmitriy Paliy. Tampere University of Technology. TICSP. 16-02-2005
% dmitriy.paliy@tut.fi

if nargin<6, derivative = 0; end;
if nargin<5, WindowType = 'Gaussian'; end;
if nargin<4, m=0; end;

xy1 = squeeze(x{1}(:)); % can be of different length
xy2 = squeeze(x{2}(:));

xz1 = squeeze(z(1,:));
xz2 = squeeze(z(2,:));
z   = squeeze(z(3,:));

sizex = length(xz1);

sizexy1 = length(xy1);
sizexy2 = length(xy2);

switch m
    case 0
        M = 1;
    case 1
        M = 3;
    case 2
        M = 6;
    case 3
        M = 10;
end;

fi_h_zero = zeros([1 M]);
fi_h_zero(derivative+1) = 1;

for i1=1:sizexy1,
    for i2=1:sizexy2,
        
        for t1=1:M,
            for t2=1:M,
                
                Fi(t1,t2) = 0;
                for j1=1:sizex,
                    Fi(t1,t2) = Fi(t1,t2) +...
                        function_Window1D(xy1(i1)-xz1(j1),h,WindowType)*...
                        function_Window1D(xy2(i2)-xz2(j1),h,WindowType)*...
                        function_Phi2D((xy1(i1)-xz1(j1))/h,(xy2(i2)-xz2(j1))/h,t1)*...
                        function_Phi2D((xy1(i1)-xz1(j1))/h,(xy2(i2)-xz2(j1))/h,t2);
                end;
                
            end;
        end;
        
        Fi = pinv(Fi);

        for j1=1:sizex,
            for t=1:M,
                fi_h(t) = function_Phi2D((xy1(i1)-xz1(j1))/h,(xy2(i2)-xz2(j1))/h,t);
            end;
            
            gh(j1) = function_Window1D(xy1(i1)-xz1(j1),h,WindowType)*...
                     function_Window1D(xy2(i2)-xz2(j1),h,WindowType)*...
                     fi_h_zero*Fi*fi_h';
        end;
        
        yh(i1,i2) = sum(sum(gh.*z));
    end;
end;