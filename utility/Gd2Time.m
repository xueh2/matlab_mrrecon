function [T1, T2] = Gd2Time(T1_0, T2_0, r1, r2, Gd)
% compute T1 and T2, given Gd conc and T1_0 and T2_0
% all times are in ms
% Gd: mmol/l
% [T1, T2] = Gd2Time(T1_0, T2_0, r1, r2, Gd)

T1_0 = T1_0/1e3;
T2_0 = T2_0/1e3;

T1 = 1 / (1/T1_0 + r1*Gd);
T1 = T1*1e3;

T2 = 1 / (1/T2_0 + r2*Gd);
T2 = T2*1e3;