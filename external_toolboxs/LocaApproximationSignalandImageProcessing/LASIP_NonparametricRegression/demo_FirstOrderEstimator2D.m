% The first-order LPA estimate (m = 1) for the 2D signal z(X_s) on the 
% IRREGULAR grid X_s is presented here. Illustration of using the function 
% yh = function_2DEstimator(...);
%
% SYNTAX
%   yh = function_2DEstimator(z,x,h,m,WindowType,derivative);
%   yh = function_2DEstimator(z,x,h,m,WindowType);
%   yh = function_2DEstimator(z,x,h,m);
%   yh = function_2DEstimator(z,x,h);
% 
% DESCRIPTION
%   The function returns the LPA estimate of a signal or its derivative. 
%   The grid of arguments X_s where observations z are made is irregular in 
%   general. The grid of arguments x where the estimate yh s done is 
%   irregular as well. The estimate is obtained by the first-order (m=1) 
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
%     m - is an order of estimator bandwidth (any real). In this example 
%         m = 1.
%
%     WindowType: 'Gaussian' (default), 'GaussianLeft', 'GaussianRight',
%                 'Rectangular', 'RectangularLeft', 'RectangularRight',
%                 'Treecube', 'Hermitian',
%                 'Exponential', 'ExponentialLeft', 'ExponentialRight',
%                 'InterpolationWindow1', 'InterpolationWindow2', 
%                 'InterpolationWindow3'
%
%     derivative - (r in the book notations) is an order of estimation 
%         derivative, which is 0 in this example.
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

clear all

% ---------------------------------------------------
% settings of the estimator
% ---------------------------------------------------
m = 1;

h = 1/8;

WindowType = 'Gaussian';

% ---------------------------------------------------
% observations
% ---------------------------------------------------
% irregular grid with function z defined on this grid
% arguments X_s
X1_s = [0.15    0.15    0.15    0.5     0.85    0.85   0.85];
X2_s = [0.3     0.5     0.7     0.5     0.3     0.5    0.7];

z(1,:) = X1_s;
z(2,:) = X2_s;
z(3,:) = [2.0     1.5    1.0      1.5     1.0     1.5     2.0]; % function values

% ---------------------------------------------------
% estimation at...
% ---------------------------------------------------
% rectangular grid for estimation
% however, it can be irregular
x1 = 0.0:1/63:1;
x2 = 0.0:1/31:1;

x{1} = x1;
x{2} = x2;

% ---------------------------------------------------
% get result
% ---------------------------------------------------
% there is NO derivative
yh = function_Estimator2D(z,x,h,m,WindowType);

% ---------------------------------------------------
% show result
% ---------------------------------------------------
version -release; % get matlab release
matlab_R=str2num(ans);

figure, mesh(x2,x1,yh),
title(['First-order approximation with ',WindowType,' window and h=',num2str(h)]);

% ---------------------------------------------------
% show observations
% ---------------------------------------------------
for i=1:length(X1_s),
    line([X2_s(i) X2_s(i)], [X1_s(i) X1_s(i)], [0 z(3,i)],...
        'Marker','o',...
        'LineStyle','--',...
        'Color','black');
end;

view(-25.50,48);

if matlab_R>=14,
    xlabel('$x_1$','interpreter','latex');
    ylabel('$x_2$','interpreter','latex');
    zlabel('$\hat{y}_h(x_1,x_2)$','interpreter','latex');
else
    xlabel('\itx_1');
    ylabel('\itx_2');
    zlabel('\ity_{\ith}(\itx_1,\itx_2)');
end;