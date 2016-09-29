function [S] = magnitude_inversion_recovery(params, TI)
% function [S]= magnitude_inversion_recovery(params, TI)
%
% function to calculate magnitude inversion recovery data (S)
% at given inversion times (TI) for specified signal amplitude and tissue
% T1.

A=params(1);
T1=params(2);

% Magnitude Inversion Recovery
S = A*abs((1 - 2*exp(-TI./T1)));