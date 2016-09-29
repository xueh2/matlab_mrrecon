% Calculate the noise covariance matrix from the noise scan

function Psi = Noise_Cov(Noise)

if size(Noise, 2) > size(Noise, 1), Noise = Noise'; end
% SN = size(Noise);
% N_std = sqrt(var(Noise));  
% 
% Noise_Scaled = Noise./(ones(SN(1), 1)*N_std);
% N_std_Scaled = N_std/min(N_std);
% 
% N_Corr = Noise_Scaled'*Noise_Scaled/SN(1); 
% 
% [V,D]=eig(N_Corr); 
% E = diag( D);
% [Cutoff, Variance_0, ks, beta,p_value, H] = KS_Cutoff_2Steps(E, SN(1) );
% for i=1:Cutoff
%     D(i,i)=1;
% end
% Psi = (N_std_Scaled'*N_std_Scaled).*(V*D*V');
% Psi = (Psi + Psi')/2;

Psi = cov(Noise);

