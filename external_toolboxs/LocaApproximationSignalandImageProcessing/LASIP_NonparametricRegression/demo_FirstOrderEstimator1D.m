% The first-order LPA estimate (m = 1) for the 1D noisy signal z(X_s) =
% y(X_s) + noise(X_s) on the IRREGULAR grid X_s is presented here. 
% Illustration of using the function [yh, yh_d1] = 
% function_FirstOrderEstimator(...).
%
% SYNTAX
%   [yh, yh_d1] = function_FirstOrderEstimator(z,x,h,WindowType);
%   [yh, yh_d1] = function_FirstOrderEstimator(z,x,h);
% 
% DESCRIPTION
%   The function returns the LPA estimate of y(X_s) from its noisy 
%   observation z(X_s). The grid of arguments X_s where observations z are 
%   made is irregular in general. The grid of arguments x where the estimate 
%   yh of y is done is irregular as well. The estimate is obtained by the 
%   first-order (m=1) approximation with a scale h.
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
% RETURNS
%   yh - is an estimated signal
% 
%   yh_d1 - is a fisrt derivative of an estimated signal
%
% REMARKS
%   The example is illustarated in Figure 2.9 of the book. For more details
%   read 2.2.5 Examples part pp. 38-40.
%
% Dmitriy Paliy. Tampere University of Technology. TICSP. 16-02-2005
% dmitriy.paliy@tut.fi

clear all

% ---------------------------------------------------
% settings of the estimator
% ---------------------------------------------------
WindowType = 'Gaussian';

h = 1/30;

% ---------------------------------------------------
% observations
% ---------------------------------------------------
X_s = [0.2139 0.3200  0.4386 0.4983 0.6435 0.7889 0.9901];

z(1,:) = X_s;
z(2,:) = [0.0185 0.82141 0.4447 0.6154 0.7919 0.9218 0.7382];

% ---------------------------------------------------
% estimation at...
% ---------------------------------------------------
x = [0.2:1/255:1];


% ---------------------------------------------------
% get result
% ---------------------------------------------------
[yh, yh_d1] = function_FirstOrderEstimator(z,x,h,WindowType);

% ---------------------------------------------------
% show result
% ---------------------------------------------------
version -release; % get matlab release
matlab_R=str2num(ans);

figure, title(['First-order approximation with ',WindowType,' window and h=',num2str(h)]),
hold on, plot(z(1,:),z(2,:),'o'), hold on, plot(x,yh,'r-'), axis tight;
if matlab_R>=14,
    xlabel('$x$','interpreter','latex');
    ylabel('$\hat{y}_h(x)$','interpreter','latex');
else
    xlabel('\itx');
    ylabel('\ity_{\ith}(\itx)');
end;