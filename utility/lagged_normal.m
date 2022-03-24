function y = lagged_normal(x,sigma, x0, tau, ampl)
% function y = lagged_normal(x,sigma, x0, tau, ampl)
%
% The lagged_normal distribution is calculated as convolution between
%   a normal distribution (x0, sigma) and a monoexponential (time const
%   tau)

%   Reference:
%      [1] Bassingthwaighte, JA. Applications of the Lagged Normal Density
%          Curve as a Model for Arterial Dilution Curves. Circulation
%          Research, Vol XVIII, April 1966, pp 398-415.



if nargin<1
    error(message('stats:lagged_normal:TooFewInputs'));
end
if nargin < 2
    sigma = 1; x0 = 0; tau = 1;
end
if nargin < 3
    sigma = 1;
end


try
    h1 = ampl*(1/(2*pi)^.5/sigma)* exp(-(1/2)*((x - x0)/sigma).^2);
    h2 = (1/tau)*exp(-x/tau);
    y = conv(h1,h2);
    y = y(1:length(h1));
catch
    error(message('stats:lagged_normal:InputSizeMismatch'));
end




