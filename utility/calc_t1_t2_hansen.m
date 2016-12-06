function [T1, T2] = calc_t1_t2_hansen(agarose_percent,NiCl2_mM, D2O_percent)
r1 = 3.750e-4 + 8.790e-6*agarose_percent + 6.683e-4*NiCl2_mM - 2.70E-04*0.01*D2O_percent;
r2 = 1.645e-4 + 7.622e-3*agarose_percent + 7.201e-4*NiCl2_mM - 4.10E-04*0.01*D2O_percent;
T1 = 1/r1;
T2 = 1/r2;
