function [yh] = function_ZeroOrderEstimator(z,x,h,WindowType)

% The zero-order LPA or Nadaray-Watson estimate (m = 0) for the 1D noisy 
% signal z(X_s) = y(X_s) + noise(X_s) on the IRREGULAR grid X_s.
%
% SYNTAX
%   yh = function_ZeroOrderEstimator(z,x,h,WindowType);
%   yh = function_ZeroOrderEstimator(z,x,h);
% 
% DESCRIPTION
%   The function returns the LPA estimate of y(X_s) from its noisy 
%   observation z(X_s). The grid of arguments X_s where observations z are 
%   made is irregular in general. The grid of arguments x where the estimate 
%   yh of y is done is irregular as well. The estimate is obtained by the 
%   zero-order (m=0) approximation with a scale h.
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
% REMARKS
%   The example is illustarated on Figure 2.6 of the book. For more details
%   read 2.2.5 Examples part pp. 34-38.
%
% Dmitriy Paliy. Tampere University of Technology. TICSP. 16-02-2005
% dmitriy.paliy@tut.fi

if nargin<4, WindowType = 'Gaussian'; end;
    
xz = squeeze(z(1,:));
z = squeeze(z(2,:));


gh = zeros(size(z));
yh = zeros(size(z));

for i=1:length(x),

    sum1 = 0;
    for s=1:length(xz),
        sum1 = sum1 + function_Window1D(x(i)-xz(s),h,WindowType);
    end;
    
    for j=1:length(xz),
        if sum1~=0,
            gh(j) = function_Window1D(x(i)-xz(j),h,WindowType)/sum1;
        else
            gh(j) = 0;
        end;
    end;
    
    yh(i) = sum(gh.*z);
end;