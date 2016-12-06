

function [T1, T2] = calc_t1_t2_captur(agarose, NiCl2)
r1 = 3.750e-4 + 8.790e-6*agarose + 6.683e-4*NiCl2;
r2 = 1.645e-4 + 7.622e-3*agarose + 7.201e-4*NiCl2;
T1 = 1/r1;
T2 = 1/r2;
