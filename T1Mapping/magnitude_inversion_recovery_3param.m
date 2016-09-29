function [S] = magnitude_inversion_recovery_3param(params, TI)
% function [S]= magnitude_inversion_recovery(params, TI)
%
% function to calculate magnitude inversion recovery data (S)
% at given inversion times (TI) for specified signal amplitude and tissue
% T1.

A=params(1);
B=params(2);
T1_star=params(3);

% Magnitude Inversion Recovery
S = abs(A - B*exp(-TI./T1_star)); % 3 parameter model