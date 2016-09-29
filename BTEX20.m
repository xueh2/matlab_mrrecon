% BTEX20 Models a tissue cylinder consisting of two regions: plasma, and
% interstitial fluid (ported from Physiome JSim model 080 by DH  9/1/12).
function BTEX20
  close all

  % Parameters
  Fp    = 1/60;   % ml/g/sec
  hVp   = 0.05;   % ml/g
  hVisf = 0.15;   % ml/g
  PSg   = .5/60;  % ml/g/sec
  Gp    = 0;      % ml/g/sec
  Gisf  = 0;      % ml/g/sec
  Dp    = 1e-5;   % cm^2/sec
  Disf  = 1e-6;   % cm^2/sec

  tmin  = 0;      % sec
  tmax  = 30;     % sec
  tdelt = 0.05;   % sec

  xmin  = 0;      % cm
  L     = 0.1;    % cm
  xdelt = L/60;   % cm

  % Spacing in x, t
  xspan = xmin:xdelt:L;
  tspan = tmin:tdelt:tmax;

  % Solver Parameters
  options = ddeset('RelTol', 1e-7, 'AbsTol', 1e-7);

  % Inflowing Concentration
  Cin = LagNormal(5, .3, 1.3, tspan);

  sol = pdepe(0, @pdefun, @icfun, @bcfun, xspan, tspan, options);
  
  figure('name', 'Cin_Cout', 'NumberTitle', 'off', 'Units', 'normalize', ...
         'Position', [0 0 .5 1])
  subplot(2, 1, 1)
  plot(tspan, Cin, 'k', tspan, sol(:, end, 1), 'r', tspan, sol(:, end, 2), 'm')
  title('Inflow, Outflow, and ISF Concentrations')
  xlabel('Time, sec')
  ylabel('Concentration mmol/ml')
  legend('Cin mM', 'Cout mM', ['Cisf(t,' num2str(L) ') mM'])

  subplot(2, 1, 2)
  plot(xspan, sol(find(tspan == 3.5), :, 1), 'k', ...
       xspan, sol(find(tspan == 4), :, 1), 'r', ...
       xspan, sol(find(tspan == 5), :, 1), 'm', ...
       xspan, sol(find(tspan == 5.5), :, 1), 'y', ...
       xspan, sol(find(tspan == 6), :, 1), 'g', ...
       xspan, sol(find(tspan == 6.5), :, 1), 'b', ...
       xspan, sol(find(tspan == 6.5), :, 2), 'b.');
  title('X-Profiles')
  xlabel('Distance along capillary, cm')
  ylabel('Concentration mmol/ml')
  legend('Cp(3.5,x) mM', 'Cp(4.0,x) mM', 'Cp(5.0,x) mM', ...
         'Cp(5.5,x) mM', 'Cp(6.0,x) mM', 'Cp(6.5,x) mM', 'Cisf(6.5,x) mM')
  
  function [c, f, s] = pdefun(x, t, u, DuDx)
    c = ones(2, 1);
    f = [Dp; Disf] .* DuDx;
    s = [-Fp * L / hVp * DuDx(1); 0] - [Gp / hVp; Gisf / hVisf] .* u + ...
        [PSg / hVp; -PSg / hVisf] .* (u(2) - u(1));
  end % pdefun

  function u = icfun(x)
    u = zeros(2, 1);
  end % icfun

  function [pl, ql, pr, qr] = bcfun(xl, ul, xr, ur, t)
    cint = interp1(tspan, Cin, t);
    pl = [-Fp * L / hVp * (ul(1) - cint); 0];
    ql = ones(2, 1);
    pr = zeros(2, 1);
    qr = ones(2, 1);
  end % bcfun
end % BTEX20

% LagNormal Lagged Normal Density Function.
% pdf = LagNormal(tSpan, tMean, RD, skewn) returns a vector of the same
% size as tSpan with the relative likelihoods that a random variable will
% take the values of a Lagged Normal Probability Density Function. The PDF
% is characterized by the parameters tMean (mean transit time), RD
% (relative dispersion) and skewn(skewness) as described in the reference.
% tSpan should be a monotonously increasing or decreasing vector. For a
% sufficiently long interval that includes the peak, the area under the
% density function approaches 1.
%
% Reference:
% Bassingthwaighte, J.B., F.H. Ackerman and E.H. Wood: Applications of the
% lagged normal density curve as a model for arterial dilution curves.
% Circ. Res. 18:398-415, 1966.
function pdf = LagNormal(tMean, RD, skewn, tspan)
  tau = RD * tMean * (skewn / 2)^(1 / 3);
  sigma = sqrt((RD * tMean) ^ 2 - tau ^ 2);
  mean = tMean - tau;
  
  options = ddeset('RelTol', 1e-7, 'AbsTol', 1e-7);
  sol =  ode23(@dCdt, tspan, 0, options);
  pdf = interp1(sol.x, sol.y, tspan);
    
  function Ret = dCdt(t, C)
    Cin = exp(-.5 * ((t - mean) / sigma) ^ 2) / (sigma * sqrt(2 * pi));
    Ret = (Cin - C) / tau;
  end % dCdt
end % LagNormal
