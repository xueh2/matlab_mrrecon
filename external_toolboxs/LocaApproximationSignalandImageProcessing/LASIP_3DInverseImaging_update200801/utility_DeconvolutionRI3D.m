% REGULARIZED INVERSE
%
% The Regularized Inverse (RI) reconstruction of slices of a 3D object at 
% specified positions and their denoising by means of the LPA-ICI technique.
% 
% Dmitriy Paliy. Tampere University of Technology.
% Updated 31-01-2008
% dmitriy.paliy@tut.fi

htimebar = timebar('RI Loop Counter','Progress'); itc=0;
totalcounter = 2*length(SReconstructed);
Precooked = 1;
% RI reconstruction
[YRIestFFT, PsiVVPsi] = function_RI3D(Z,V,RIreg2,SObserved,SReconstructed);

h1 = figure('name','RI Estimates');

for sorig=1:length(SReconstructed),
    itc = itc+1;  timebar(htimebar,itc/totalcounter);

    sigma = sigmanoise_arr(sorig);
    EY{sorig} = real(ifft2(squeeze(YRIestFFT(sorig,:,:))));

    if DoFiltering==1,
        stdh = squeeze(PsiVVPsi(sorig,sorig,1:xN,1:xN));
        sigma = sigmanoise_arr(sorig);

        % LPA-ICI denoising for each slice
        [YRIest(sorig,:,:), SigmaICI] = function_LPAICIDenoising(fft2(EY{sorig}),stdh,sigma,GammaParameterRI,Precooked);
        EY{sorig} = squeeze(YRIest(sorig,:,:));
        EY{sorig} = min(EY{sorig},1); EY{sorig} = max(EY{sorig},0);
        YRIest_All{sorig}(:,:) = EY{sorig};
        SigmaICI_All{sorig}(:,:) = SigmaICI;
    else
        EY{sorig} = min(EY{sorig},1); EY{sorig} = max(EY{sorig},0);
        SigmaICI_All{sorig}(:,:) = squeeze(PsiVVPsi(sorig,sorig,:,:));
        YRIest_All{sorig}(:,:) = EY{sorig};
    end;

    itc = itc+1;  timebar(htimebar,itc/totalcounter);

    last_errors = function_Errors(Y{sorig},EY{sorig},Z{1});
    
    figure(h1),
    subplot(1,length(SReconstructed),sorig), imshow(EY{sorig},[]),
    title(['RI Est. RMSE = ',num2str(last_errors(5)),...
        '; PSNR = ',num2str(last_errors(3))]);

    drawnow;
end;

close(htimebar);