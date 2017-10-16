function R = AATH_s(Fp, Vp, PS, Ve, s)
% compute AATH model
% R = AATH_s(Fp, Vp, PS, Ve, s)
% R = ( 1-exp(-s*Tc) )/s + E*exp(-s*Tc)/(kep+s)
% Tc = Vp/Fp
% kep = EFp/Ve
% E = 1-exp(-PS/Fp)

Tc = Vp/Fp;
E = 1-exp(-PS/Fp);
kep = E*Fp/Ve;                                                                             

ind = find(abs(s)==0);

% R = (1-exp(-s(:)*Tc))./s(:) + E*exp(-Tc*s(:))./(s(:)+kep);
% R(ind(:)) = Tc + E/kep;

R = (1-exp(-s(:)*Tc))./s(:) + E*exp(-Tc*s(:))./(s(:)+kep);
R(ind(:)) = Tc + E/kep;