function [S, J] = phase_sensitive_inversion_recovery_3param(params, TI)

A=params(1);
B=params(2);
T1_star=params(3);

% Phase Sensitive Inversion Recovery (PSIR)
S = A - B*exp(-TI./T1_star); % 3 parameter model

if nargout == 2
   % partial derivatives (Jacobian) with respect to params(i) 
   J = [(1 - B*exp(-TI./T1_star)); -B*exp(-TI/T1_star);...
         -B/(T1_star^2).*exp(-TI/T1_star)]';
end