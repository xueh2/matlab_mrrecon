function [F, tau, k, td, yr, r] = PerformDeconvolution_Fermi(cin, y, F0, tau0, k0, td0)
% perform the fermi deconvolution

%% set up common parameters
N = numel(cin);

if numel(y)~=N
    error('cin and y have different length');
end

if nargin < 3
    F0 = 0.015;
end

if nargin < 4
    tau0 = 1.0;
end

if nargin < 5
    k0 = 0.1;
end

if nargin < 6
    td0 = 0;
end

%% compute the fermi deconvolution

A = zeros(N, N);
for i=1:N
    for j=i:-1:1
        A(i,j) = cin(i-j+1);
    end
end

t = [0:N-1]';

FAll = F0;
tauAll = 0.0001:0.05:2;
kAll = -1:0.02:2;

cost = zeros(numel(FAll), numel(tauAll), numel(kAll));

for f=1:numel(FAll)
    for m=1:numel(tauAll)
        for n=1:numel(kAll)
            cost(f, m, n) = fermi([FAll(f), tauAll(m) kAll(n)]);
        end
    end
end

[mc, ind] = min(cost(:));
[I, J, S] = ind2sub(size(cost), ind);

F0 = FAll(I);
tau0 = tauAll(J);
k0 = kAll(S);

[X,FVAL,EXITFLAG,OUTPUT] = fminsearch(@fermi, [F0 tau0 k0]);

F = X(1);
tau = X(2);
k = X(3);
td = 0;

r = F ./ (1 + exp( k * ( tau-t ) ) );

yr = A*r;

    function v = fermi(par)
        % compute impulse                      
        r = par(1) ./ (1 + exp( par(3) * ( par(2)-t ) ) );
        
        diff = reshape(y - A*r, [N 1]);
        v = diff'*diff;
    end

end


