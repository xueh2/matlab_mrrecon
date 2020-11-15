% Dmitriy Paliy. Tampere University of Technology. 
% Updated 31-01-2008
% dmitriy.paliy@tut.fi

for i=1:length(SOriginal), Y{i} = Fxyz(:,:,SOriginal(i)); end;

sigmaz = 4;

for sobs=1:length(SObserved),
    for sorig=1:length(SOriginal),
        ALPHA(sobs,sorig) = exp(-(abs(SOriginal(sorig)-SObserved(sobs)))/(2*sigmaz^2));
    end;
    ALPHA(sobs,:) = ALPHA(sobs,:)./sum(ALPHA(sobs,:));
end;

k_coef = 1;
for sobs=1:length(SObserved),

    for sorig=1:length(SOriginal),
        sigmagh = abs(SOriginal(sorig)-SObserved(sobs))*k_coef;
        if sigmagh==0, sigmagh = eps; end;
        gh = fspecial('gaussian',11,sigmagh);
        gh = ALPHA(sobs,sorig)*gh;
        
        % gh = fspecial('average',9);
        
        [ghx,ghy] = size(gh);


        ghx2 = (ghx-1)/2;
        ghy2 = (ghy-1)/2;
        Convh_mat = zeros(xN,xN);
        Convh_mat(1:ghx,1:ghy) = gh;
        D1 = diag(exp(j*2*pi*(ghx2)*[0:xN-1]/(xN)));
        D2 = diag(exp(j*2*pi*(ghy2)*[0:xN-1]/(xN)));
        V{sobs,sorig} = D1*fft2(Convh_mat)*D2;
        z = real(ifft2(fft2(Y{sorig}).*V{sobs,sorig}));

        Ztmp(sorig,1:size(z,1),1:size(z,2)) = z(:,:);
    end;

    org_sigma = sqrt(norm(Ztmp(:)-mean(Ztmp(:)),2)^2 /(xN*yN*10^(BSNRdB/10))); % STD of Noise in Blurred Image
    sigmanoise_arr(sobs) = org_sigma;
    randn('state',2055615788+sobs);
    noise = sigmanoise_arr(sobs)*randn(xN,xN);
    cs = squeeze(sum(Ztmp,1)) + noise;

    Z{sobs} = cs;
end;

clear ALPHA;
for sobs=1:length(SObserved),
    for sorig=1:length(SOriginal),
        ALPHA(sobs,sorig) = exp(-(abs(SOriginal(sorig)-SObserved(sobs)))/(2*sigmaz^2));
    end;
    ALPHA(sobs,:) = ALPHA(sobs,:)./sum(ALPHA(sobs,:));
end;

for sobs=1:length(SObserved),
    for sorig=1:length(SReconstructed),
        sigmagh = abs(SReconstructed(sorig)-SObserved(sobs))*k_coef;
        if sigmagh==0, sigmagh = eps; end;
        gh = fspecial('gaussian',11,sigmagh);
        gh = gh./sum(gh(:));
        gh = ALPHA(sobs,sorig).*gh;
        
        % gh = fspecial('average',9);

        [ghx,ghy] = size(gh);

        ghx2 = (ghx-1)/2;
        ghy2 = (ghy-1)/2;
        Convh_mat = zeros(xN,xN);
        Convh_mat(1:ghx,1:ghy) = gh;
        D1 = diag(exp(j*2*pi*(ghx2)*[0:xN-1]/(xN)));
        D2 = diag(exp(j*2*pi*(ghy2)*[0:xN-1]/(xN)));
        VR{sobs,sorig} = D1*fft2(Convh_mat)*D2;
    end;
end;

% save InputData Z V VR Y;
