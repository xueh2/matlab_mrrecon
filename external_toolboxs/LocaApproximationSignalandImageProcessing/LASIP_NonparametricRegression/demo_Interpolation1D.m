% Demo of using 1D Interpolation Estimator on irregular grid of orders m = 0 
% amd m = 1 with interpolation wondows.
%
% SYNTAX
%   yh = function_ZeroOrderEstimator(z,x,h,WindowType);
%   yh = function_ZeroOrderEstimator(z,x,h);
%
% DESCRIPTION
%   The function returns the LPA estimate of y(X_s) from its observation 
%   z(X_s) using interpolation window. The grid of arguments X_s is 
%   irregular in general. The grid of arguments x where the interpolation 
%   estimate yh of y is done is irregular as well. The estimate is obtained 
%   by the zero-order (m=0) and first-order (m=1) approximation with scales h.
%
%   Here:
%     z - is a 2xn array consists of:
%       a) z(1,1:n) is a grid of arguments X_s where function z is defined. 
%          It has length n and can be regular or irregular;
%       b) z(2,1:n) is an array of function z values at points X_s.
% 
%     x - is an array of ARGUMENTS where the LPA estimation is done and, in 
%         general, is different from X_s and also irregular.
% 
%     h - estimator bandwidth (any real);
%
%     WindowType: 'Gaussian' (default), 'GaussianLeft', 'GaussianRight',
%                 'Rectangular', 'RectangularLeft', 'RectangularRight',
%                 'Treecube', 'Hermitian',
%                 'Exponential', 'ExponentialLeft', 'ExponentialRight',
%                 'InterpolationWindow1', 'InterpolationWindow2', 
%                 'InterpolationWindow3'
%
%   In this example WindowType = 'InterpolationWindow1' is used.
%
% REMARKS
%   The example is illustarated on Figure 2.13 of the book. For more details
%   read part 2.5 'Nonparametric interpolation' pp. 47-49.
%
% Dmitriy Paliy. Tampere University of Technology. TICSP. 16-02-2005
% dmitriy.paliy@tut.fi

clear all

% ---------------------------------------------------
% settings of the estimator
% ---------------------------------------------------
WindowType = 'InterpolationWindow1';

h = 1/120;

% ---------------------------------------------------
% observations
% ---------------------------------------------------
X_s = [0.21   0.32    0.4983   0.6435  0.99];

z(1,:) = X_s;
z(2,:) = [0.7621 0.45647 0.018504 0.82141 0.4447]; 

% ---------------------------------------------------
% estimation at regular grid
% ---------------------------------------------------
x = [0:(1.2/255):1.2];

% ---------------------------------------------------
% get result for zero-order
% ---------------------------------------------------
yh0 = function_ZeroOrderEstimator(z,x,h,WindowType);

% ---------------------------------------------------
% get result for first-order
% ---------------------------------------------------
yh1 = function_FirstOrderEstimator(z,x,h,WindowType);

% ---------------------------------------------------
% show results
% ---------------------------------------------------
version -release; % get matlab release
matlab_R=str2num(ans);

figure, 
subplot(1,2,1), plot(X_s,z(2,:),'o'), hold on, plot(x,yh0,'r-'), axis tight,
title(['Zero-order interpolation with ',WindowType,' window and h=',num2str(h)]);
if matlab_R>=14,
    xlabel('$x$','interpreter','latex');
    ylabel('$\hat{y}_h(x)$','interpreter','latex');
else
    xlabel('\itx');
    ylabel('\ity_{\ith}(\itx)');
end;

subplot(1,2,2), plot(X_s,z(2,:),'o'), hold on, plot(x,yh1,'r-'), axis tight,
title(['First-order interpolation with ',WindowType,' window and h=',num2str(h)]);
if matlab_R>=14,
    xlabel('$x$','interpreter','latex');
    ylabel('$\hat{y}_h(x)$','interpreter','latex');
else
    xlabel('\itx');
    ylabel('\ity_{\ith}(\itx)');
end;