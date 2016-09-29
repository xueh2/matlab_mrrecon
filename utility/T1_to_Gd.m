function Gd = T1_to_Gd(T1_0, T1, r1)
% Given T1_0 and T1, compute Gd
% all times are in ms
% Gd: mmol/l
% Gd = T1_to_Gd(T1_0, T1, r1)

T1_0 = T1_0/1e3;
T1 = T1/1e3;

Gd = (1/T1 - 1/T1_0 ) / r1;