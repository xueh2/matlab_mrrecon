function [Ki, tau, k, Vb, E, yr, r] = PerformDeconvolution_TwoCompFermiExtraction(cin, y, Ki0, tau0, k0, Vb0, E0)
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

if nargin < 7
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
tauAll = 0.0001:0.05:2;
kAll = -1:0.02:2;
VbAll = 0:0.02:0.2;
EAll = 0.3:0.05:0.75;

cost = zeros(numel(KiAll), numel(tauAll), numel(kAll), numel(VbAll), numel(EAll));

for f=1:numel(KiAll)
    for m=1:numel(tauAll)
        for k=1:numel(kAll)
            for n=1:numel(VbAll)
                for e=1:numel(EAll)
                    cost(f, m, k, n, e) = two_comp_fermi_extraction([KiAll(f) tauAll(m) kAll(k) VbAll(n) EAll(e)]);
                end
            end
        end
    end
end

[mc, ind] = min(cost(:));
[I, J, S, P, Q] = ind2sub(size(cost), ind);

Ki0 = KiAll(I);
tau0 = tauAll(J);
k0 = kAll(S);
% Vb0 = VbAll(P);
% E0 = EAll(Q);
Vb0 = 0.1;
E0 = 0.5;

% [X,FVAL,EXITFLAG,OUTPUT] = fminsearch(@two_comp_fermi_extraction, [Ki0 tau0 k0 Vb0 E0]);

options = optimoptions('fmincon','Algorithm','interior-point');
[X,FVAL,EXITFLAG] = fmincon(@two_comp_fermi_extraction, [Ki0 tau0 k0 Vb0 E0], [], [], [], [], [1e-3 -1 -1 0.01 0.1], [100 100 100 0.6 1], [], options);

Ki = X(1);
tau = X(2);
k = X(3);
Vb = X(4);
E  = X(5);

r1 = Ki ./ (1 + exp( k * ( tau-t ) ) );
r2 = Vb * (1-E) .* cin;
k2 = max(r1(:)) /  sum(r1(:)); % reciprocal of mean transit time
r3 = Vb * E * k2/Ki .* r1;

r = r1;

yr = A*(r1+r3)+r2;

    function v = two_comp_fermi_extraction(par)
        % compute impulse
       
        r1 = par(1) ./ (1 + exp( par(3) * ( par(2)-t ) ) );
        r2 = par(4) * (1-par(5)) .* cin;
        k2 = max(r1(:)) /  sum(r1(:)); % reciprocal of mean transit time
        r3 = par(4) * par(5) * k2/par(1) .* r1;
       
        diff = reshape(y - A*(r1+r3)+r2, [N 1]);
        v = diff'*diff;
    end

end


