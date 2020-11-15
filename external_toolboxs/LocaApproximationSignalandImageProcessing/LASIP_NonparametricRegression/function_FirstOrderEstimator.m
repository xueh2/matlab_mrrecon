function [yh, yh1] = FirstOrderEstimator(z,x,h,WindowType)

% The first order LPA estimate (m = 1) for the 1D noisy signal z(X_s) =
% y(X_s) + noise(X_s) on the IRREGULAR grid X_s.
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
%   first order (m=1) approximation with a scale h.
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

if nargin<4, WindowType = 'Gaussian'; end;

xy = x;

xz = z(1,:);
z = z(2,:);

sizex = length(xz);

sizey = length(xy);

gh = zeros(size(xz));
gh1 = zeros(size(xz));
yh = zeros(size(xy));
yh1 = zeros(size(xy));
Sl = zeros(1,3);

for i=1:sizey,

    for j=0:2,
        Sl(j+1) = 0;
        for s=1:sizex,
            Sl(j+1) = Sl(j+1) +...
                function_Window1D(xy(i)-xz(s),h,WindowType)*...
                (xy(i)-xz(s))^j;
        end;
    end;

    for j1=1:sizex,

        devident = Sl(2+1)-Sl(1+1)*(xy(i)-xz(j1));
        devisor = Sl(2+1)*Sl(0+1)-Sl(1+1)^2;
        if devisor~=0,
            gh(j1) = function_Window1D(xy(i)-xz(j1),h,WindowType)*...
                devident/devisor;
        else
            gh(j1) = 0;
        end;

        devident1 = -Sl(1+1)+Sl(0+1)*(xy(i)-xz(j1));
        devisor1 = Sl(2+1)*Sl(0+1)-Sl(1+1)^2;
        if devisor1~=0,
            gh1(j1) = -function_Window1D(xy(i)-xz(j1),h,WindowType)*...
                devident1/devisor1;
        else
            gh1(j1) = 0;
        end;        
    end;

    yh(i) = sum(gh.*z);

    yh1(i) = sum(gh1.*z);
end;