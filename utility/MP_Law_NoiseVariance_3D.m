
function [Cutoff, Variance, SNR, gfactor] = MP_Law_NoiseVariance_3D(data)
% This is a function to estimate the noise variance in 3D image using the MP Law.
% [Cutoff, Variance] = MP_Law_NoiseVariance_3D(data);
% data: original data matrix [Nfe Npe frame]
     
s = size(data);
[V, D] = KL_Eigenvalue(data);
E = diag(D); %figure(4), hist(E, 24)
[Cutoff, Variance, ks, beta, p_value, H] = KS_Cutoff_2Steps(E, s(1)*s(2));

SNR = sqrt(max(E)/s(3))/sqrt(Variance);

%[a,V,D, eigImg] = KL_Eigenimage_NoMean(data);
[eigImg, V, D] = KL_Eigenimage(data);

gfactor = std(eigImg(:,:,1:Cutoff),0,3);