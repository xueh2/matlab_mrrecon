function [Ki, k2, Vb, yr, r] = PerformDeconvolution_TwoCompExp(cin, y, withVb, Ki0, k20, Vb0)
% perform the fermi deconvolution

%% set up common parameters
N = numel(cin);

if numel(y)~=N
    error('cin and y have different length');
end

if nargin < 4
    Ki0 = 0.015;
end

if nargin < 5
    k20 = Ki0 * (1-0.45)/0.25;
end

if nargin < 6
    Vb0 = 0.1;
end

%% compute the fermi deconvolution

A = zeros(N, N);
for i=1:N
    for j=i:-1:1
        A(i,j) = cin(i-j+1);
    end
end

t = [0:N-1]';

KiAll = Ki0;
k2All = Ki0 * (1-0.45)./[0.01:0.01:0.3];
VbAll = 0:0.01:0.2;

cost = zeros(numel(KiAll), numel(k2All), numel(VbAll)) + 1e8;

for f=1:numel(KiAll)
    for m=1:numel(k2All)
        for n=1:numel(VbAll)
            cost(f, m, n) = two_comp_exp([KiAll(f), k2All(m) VbAll(n)]);
        end
    end
end

[mc, ind] = min(cost(:));
[I, J, S] = ind2sub(size(cost), ind);

Ki0 = KiAll(I);
k20 = k2All(J);
Vb0 = VbAll(S);

options = optimset('MaxFunEvals', 1000);
tolx = 1e-7;
tolf = 1e-7;
maxiter = 1000;
maxfun = 600;
prnt = 0;

if(withVb)
    [X,FVAL,EXITFLAG,OUTPUT] = fminsearch_mr(@two_comp_exp, [Ki0 k20 Vb0], tolx, tolf, maxiter, maxfun, prnt);
    %  [X,FVAL,EXITFLAG,OUTPUT] = fminsearch(@two_comp_exp, [Ki0 k20 Vb0]);

    % options = optimoptions('fmincon','Algorithm','interior-point');
    % [X,FVAL,EXITFLAG] = fmincon(@two_comp_exp, [Ki0 k20 Vb0], [], [], [], [], [1e-3 0 0], [1 1 1], [], options);

    Ki = X(1);
    k2 = X(2);
    Vb = X(3);

    % r =  (1-Vb) * Ki .* exp(-k2.*t) + Vb;
    r =  Ki .* exp(-k2.*t) + Vb;
else
    [X,FVAL,EXITFLAG,OUTPUT] = fminsearch_mr(@two_comp_exp2, [Ki0 k20], tolx, tolf, maxiter, maxfun, prnt);
    
    Ki = X(1);
    k2 = X(2);
    Vb = 0;
    r =  Ki .* exp(-k2.*t);
end

yr = A*r;

    function v = two_comp_exp(par)
        % compute impulse                      
        % r = (1-par(3)) * par(1) .* exp(-par(2).*t) + par(3);
        % r = par(1) .* exp(-par(2).*t) + par(3);
        r = abs(par(1)) .* exp(- abs(par(2)).*t) + abs(par(3));
        
        diff = reshape(y - A*r, [N 1]);
        v = diff'*diff;
    end

    function v = two_comp_exp2(par)
        % compute impulse                      
        % r = (1-par(3)) * par(1) .* exp(-par(2).*t) + par(3);
        % r = par(1) .* exp(-par(2).*t) + par(3);
        r = abs(par(1)) .* exp(- abs(par(2)).*t);
        
        diff = reshape(y - A*r, [N 1]);
        v = diff'*diff;
    end

end


