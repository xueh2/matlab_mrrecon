function T1 = Gd2T1(Gd, T1_0, r1);
% function T1 = Gd2T1(Gd, T1_0, r1);
%
% compute T1 (sec) from [Gd] (mM), T1_0 (sec), and r1

R1 = (1/T1_0) + r1*Gd;
T1 = 1./R1;