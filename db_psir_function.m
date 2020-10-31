function M = db_psir_function(T1, T2, TD1, TD2, TE);
% function M = db_psir_function(T1, T2, TD1, TD2, TE);
%
% M = db_psir_function(T1, T2, TD1, 0, 0); % becomes normal PSIR

M1 = (1 - 2*exp(-TD1./T1)).*exp(-TE./T2);
M = 1 - (1 - M1).*exp(-TD2./T1);
            
