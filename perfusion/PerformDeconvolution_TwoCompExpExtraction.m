function [Ki, k2, Vb, E, yr, r] = PerformDeconvolution_TwoCompExpExtraction(cin, y, Ki0, k20, Vb0, E0)
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
    k20 = Ki0 * (1-0.45)/0.25;
end

if nargin < 5
    Vb0 = 0.1;
end

if nargin < 5
    E0 = 0.6;
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
EAll = 0.5:0.02:1;

cost = zeros(numel(KiAll), numel(k2All), numel(VbAll), numel(EAll));

for f=1:numel(KiAll)
    for m=1:numel(k2All)
        for n=1:numel(VbAll)
            for e=1:numel(EAll)
                cost(f, m, n, e) = two_comp_exp_extraction([KiAll(f), k2All(m) VbAll(n) EAll(e)]);
            end
        end
    end
end

[mc, ind] = min(cost(:));
[I, J, S, P] = ind2sub(size(cost), ind);

Ki0 = KiAll(I);
k20 = k2All(J);
% Vb0 = VbAll(S);
% E0 = EAll(P);
Vb0 = 0.1;
E0 = 0.5;

% [X,FVAL,EXITFLAG,OUTPUT] = fminsearch(@two_comp_exp_extraction, [Ki0 k20 Vb0 E0]);

options = optimoptions('fmincon','Algorithm','interior-point');
[X,FVAL,EXITFLAG] = fmincon(@two_comp_exp_extraction, [Ki0 k20 Vb0 E0], [], [], [], [], [0 0 0.01 0.1], [100 100 0.6 1], [], options);

Ki = X(1);
k2 = X(2);
Vb = X(3);
E  = X(4);

r1 = Ki .* exp(-k2.*t);
r2 = Vb * (1-E) .* cin;
r3 = Vb * E * k2/Ki .* r1;

r = r1;

yr = A*(r1+r3)+r2;

    function v = two_comp_exp_extraction(par)
        % compute impulse
        if(par(1)<1e-6)
            par(1) = 1e-6;
        end
        
        r1 = par(1) .* exp(-par(2).*t);
        r2 = par(3) * (1-par(4)) .* cin;
        r3 = par(3) * par(4) * par(2)/par(1) .* r1;
       
        diff = reshape(y - A*(r1+r3)+r2, [N 1]);
        v = diff'*diff;
    end

end


