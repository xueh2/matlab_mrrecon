function [S] = saturation_recovery(params, TI)
% function [S]= saturation_recovery(params, TI)
%
% function to calculate magnitude saturation recovery data (S)
% at given inversion times (TI) for specified signal amplitude and tissue
% T1.

A=params(1);
T1=params(2);

% saturation recovery
S = A*abs(1 - exp(-TI./T1));