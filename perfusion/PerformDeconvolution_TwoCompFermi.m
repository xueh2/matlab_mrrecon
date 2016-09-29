function [Ki, tau, k, Vb, yr, r] = PerformDeconvolution_TwoCompFermi(cin, y, Ki0, tau0, k0, Vb0)
% perform the fermi deconvolution

%% set up common parameters
N = numel(cin);

if numel(y)~=N
    error('cin and y have different length');
end

if nargin < 3
    Ki0 = 0.015;
end

if nargin < 4
    tau0 = 1.0;
end

if nargin < 5
    k0 = 0.1;
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
tauAll = 0.0001:0.1:2;
kAll = -0.1:0.04:2;
VbAll = 0:0.02:0.2;

cost = zeros(numel(KiAll), numel(tauAll), numel(kAll), numel(VbAll));

for f=1:numel(KiAll)
    for m=1:numel(tauAll)
        for n=1:numel(kAll)
            for b=1:numel(VbAll)
                cost(f, m, n, b) = two_comp_fermi([KiAll(f), tauAll(m), kAll(n) VbAll(b)]);
            end
        end
    end
end

[mc, ind] = min(cost(:));
[I, J, S, P] = ind2sub(size(cost), ind);

Ki0 = KiAll(I);
tau0 = tauAll(J);
k0 = kAll(S);
Vb0 = VbAll(P);

[X,FVAL,EXITFLAG,OUTPUT] = fminsearch(@two_comp_fermi, [Ki0 tau0 k0 Vb0]);

% options = optimoptions('fmincon','Algorithm','interior-point', 'Display', 'off');
% [X,FVAL,EXITFLAG] = fmincon(@two_comp_fermi, [Ki0 tau0 k0 Vb0], [], [], [], [], [1e-3 -1 -1 0], [100 100 1 1], [], options);

Ki = X(1);
tau = X(2);
k = X(3);
Vb = X(4);

r = (1-Vb)* Ki ./ (1 + exp( k * ( tau-t ) ) ) + Vb;

yr = A*r;

    function v = two_comp_fermi(par)
        % compute impulse                      
        r = ( 1-par(4) ) * par(1) ./ (1 + exp( par(3) * ( par(2)-t ) ) ) + par(4);
        
        diff = reshape(y - A*r, [N 1]);
        v = diff'*diff;
    end

end


