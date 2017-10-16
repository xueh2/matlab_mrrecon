function R = AATH_t(Fp, Vp, PS, Ve, t)
% compute AATH model
% R = AATH_t(Fp, Vp, PS, Ve, t)
% R = step(Tc-t) + E*step(t-Tc)*exp( -(t-Tc)*kep )
% Tc = Vp/Fp
% kep = EFp/Ve
% E = 1-exp(-PS/Fp)

Tc = Vp/Fp;
E = 1-exp(-PS/Fp);
kep = E*Fp/Ve;                                                                             

% th1 = heaviside(Tc-t);
% th2 = heaviside(t-Tc);

% R = th1+E*th2.*exp(-(t-Tc)*kep);
R = heaviside(Tc-t)+E*heaviside(t-Tc).*exp(-(t-Tc)*kep);

% for ii=1:numel(t)
%     
%     if(t(ii)<Tc)
%         R(ii) = 1;
%     else
%         R(ii) = E*exp( -(t(ii)-Tc) * kep );
%     end
% end