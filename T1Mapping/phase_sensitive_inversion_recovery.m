function [S, J] = phase_sensitive_inversion_recovery(params, TI)
% function [S, J]= inversion_recovery(params, TI)
%
% function to calculate (phase sensitive) inversion recovery data (S)
% at given inversion times (TI) for specified signal amplitude and tissue
% T1.

A=params(1);
T1=params(2);

% Phase Sensitive Inversion Recovery (PSIR)
S = A*(1 - 2*exp(-TI./T1));

if nargout ==2
   % partial derivatives (Jacobian) with respect to paramss(i) 
   J = [(1 - 2*exp(-TI./T1)); ...
         -2*A*TI.^2/(T1^3).*exp(-TI/T1)]';
end