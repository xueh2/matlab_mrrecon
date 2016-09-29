function [F, PS, Vp, Visf, Q_e, FVAL] = PerformPDEPerfusion_BTEX20(cin, y, F0, PS0, tdelt, hVp, hVisf, Q_e, Fp, Vp, PS, Visf)
% perform the fermi deconvolution
% [F, PS, Vp, Visf, Q_e, FVAL] = PerformPDEPerfusion_BTEX20(cin, y, F0, PS0, tdelt, hVp, hVisf, Q_e, Fp, Vp, PS, Visf)

%% set up common parameters
N = numel(cin);

if numel(y)~=N
    error('cin and y have different length');
end

if nargin < 4
    tdelt = 0.5
end

if nargin < 5
    hVp = 0.05;
end

if nargin < 6
    hVisf = 0.25;
end

R = 60/tdelt;

cost = zeros(numel(Fp), numel(PS), numel(Vp), numel(Visf)) + 1e8;

for f=1:numel(Fp)
    for p=1:numel(PS)
        for v=1:numel(Vp)
            for vi=1:numel(Visf)
                qq = Q_e(:, f, v, p, vi);
                qq = qq(1:numel(y));
                cost(f, p, v, vi) = norm(qq(:) - y(:));
            end
        end
    end
end

[mc, ind] = min(cost(:));
[f, p, v, vi] = ind2sub(size(cost), ind);

F0 = Fp(f);
PS0 = PS(p);
hVp = Vp(v);
hVisf = Visf(vi);

options = optimset('MaxFunEvals', 1000);
tolx = 1e-4;
tolf = 1e-4;
maxiter = 100;
maxfun = 100;
prnt = 0;

Gp = 0;
Gisf = 0;
Dp = 1e-5;
Disf = 1e-6;
xmin  = 0;      % cm
L     = 0.1;    % cm
xdelt = L/30;   % cm
xspan = xmin:xdelt:L;

% h = figure;
% hold on
% plot(cin);
% plot(y);
        
[X,FVAL,EXITFLAG,OUTPUT] = fminsearch_mr(@BTEX20_cost, [F0 PS0 hVp hVisf], tolx, tolf, maxiter, maxfun, prnt);

F = X(1);
PS = X(2);
Vp = X(3);
Visf = X(4);

% call the BTEX20 model

    function v = two_comp_exp2(par)
        % compute impulse                      
        % r = (1-par(3)) * par(1) .* exp(-par(2).*t) + par(3);
        % r = par(1) .* exp(-par(2).*t) + par(3);
        r = abs(par(1)) .* exp(- abs(par(2)).*t);
        
        diff = reshape(y - A*r, [N 1]);
        v = diff'*diff;
    end

    function v = BTEX20_cost(par)
        % compute BTEX cost function    
        if(par(2)>3)
            par(2) = 3;
        end
        
        if(par(3)>0.1)
            par(3) = 0.1;
        end

        if(par(4)>0.5)
            par(4) = 0.5;
        end
        
        % [C_e, C_p] = BTEX20_model(cin, (numel(cin)-1)*tdelt, tdelt, par(1)/60, par(2)/60, hVp, hVisf);
        [sol_m, C_e_m, C_p, Q_e_m] = Matlab_gt_BTEX20_model( double(cin(:)), [0:(numel(cin)-1)]*tdelt, xspan, par(1), par(3), par(2), par(4), Gp, Gisf, Dp, Disf);
        
        % compute Q_e
        Q_e = zeros(numel(cin), 1);
        for tt=1:numel(cin)
            Q_e(tt) = sum(par(1)*(cin(1:tt) - C_p(1:tt)) * tdelt) / 60;
        end
            
        diff = y(:) - Q_e(:);
%         diff = y(1:valley) - Q_e(1:valley);
        v = diff'*diff;
        
%         figure(h)       
%         plot(Q_e, 'r');
    end
end


