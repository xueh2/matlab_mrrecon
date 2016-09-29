function covar = computeT1PSIRError(A, B, T1Star, TI, sigma)
% compute the T1 PSIR fitting error

if ( A < 0 )
    covar = zeros(3,3);
    return;
end

T1 = (B/A - 1) * T1Star;

N = numel(TI);

deriv = zeros(N, 3); % store derivative for every TI for A, B, T1

for n=1:N
    [dydA, dydB, dydT1] = computeDeriv(TI(n), A, B, T1);
    deriv(n,1) = dydA;
    deriv(n,2) = dydB;
    deriv(n,3) = dydT1;    
end

covar = zeros(3,3);
   
for k=1:3
    for l=1:3
        for n=1:N            
            covar(k,l) = covar(k,l) + 1.0/(sigma(n)*sigma(n)) * deriv(n,k) * deriv(n,l);
        end   
    end
end

if ( det(covar) < 1e-12 )
    covar = zeros(3,3);
    return;
end

covar = inv(covar);

function [dydA, dydB, dydT1] = computeDeriv(TI, A, B, T1)

ep = exp(-TI*(B/A-1)/T1);

dydA = 1-B*ep*TI*B/(T1*A*A);
dydB = -ep + B*ep*TI/(T1*A);
dydT1 = -B*ep*TI*(B/A-1)/(T1*T1);
