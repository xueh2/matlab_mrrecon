% ================================
% == REGULARIZED WIENER INVERSE ==
% ================================
%
% The Regularized Wiener Inverse (RWI) reconstruction of slices at 
% specified positions of a 3D object and their denoising by means of the 
% LPA-ICI technique.
%
% Dmitriy Paliy. Tampere University of Technology. 
% Updated 31-01-2008
% dmitriy.paliy@tut.fi

htimebar = timebar('RWI Loop Counter','Progress'); itc=0;
totalcounter = 2*length(SReconstructed);
Precooked = 2;
% RWI reconstruction
[YRWestFFT, PsiVVPsiRW] = function_RWI3D(Z,V,EY,sigmanoise_arr,RWreg2,SObserved,SReconstructed);

h1 = figure('name','RWI Estimates');

for sorig=1:length(SReconstructed),
    itc = itc+1;  timebar(htimebar,itc/totalcounter);

    sigma = sigmanoise_arr(sorig);
    EY_RW{sorig} = real(ifft2(squeeze(YRWestFFT(sorig,:,:))));

    if DoFiltering==1,
        stdh = squeeze(PsiVVPsiRW(sorig,sorig,1:xN,1:xN));
        sigma = sigmanoise_arr(sorig);

        % LPA-ICI denoising
        [YRWest(sorig,:,:), SigmaICI] = function_LPAICIDenoising(fft2(EY_RW{sorig}),stdh,sigma,GammaParameterRW,Precooked);
        EY_RW{sorig} = squeeze(YRWest(sorig,:,:));
        EY_RW{sorig} = min(EY_RW{sorig},1); EY_RW{sorig} = max(EY_RW{sorig},0);
        YRWest_All_RW{sorig}(:,:) = EY_RW{sorig};
        SigmaICI_All_RW{sorig}(:,:) = SigmaICI;
    else
        EY_RW{sorig} = min(EY_RW{sorig},1); EY_RW{sorig} = max(EY_RW{sorig},0);
        SigmaICI_All_RW{sorig}(:,:) = squeeze(PsiVVPsiRW(sorig,sorig,:,:));
        YRWest_All_RW{sorig}(:,:) = EY_RW{sorig};
    end;

    last_errors = function_Errors(Y{sorig},EY_RW{sorig},Z{1});
    
    figure(h1),
    subplot(1,length(SReconstructed),sorig), imshow(EY_RW{sorig},[]),
    title(['RWI Est. RMSE = ',num2str(last_errors(5)),...
        '; PSNR = ',num2str(last_errors(3))]);

    itc = itc+1;  timebar(htimebar,itc/totalcounter);
end;

close(htimebar);