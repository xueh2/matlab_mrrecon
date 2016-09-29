
function kspace = performDownSampleFE(kspace_in)
% perform the downsample for FE direction
% in siemens system, there is always 2 times oversampling along the FE

s = size(kspace_in);

Nfe = size(kspace_in, 1);

rfft = 1/sqrt(Nfe/2);
rifft = sqrt(Nfe);

if ( length(s) == 2 )       
    kspaceTemp = fftshift(kspace_in,1);
    imTemp = rifft*ifftshift(ifft(kspaceTemp, [], 1), 1);

    lenCut = Nfe / 4;    
    imTemp2 = imTemp(lenCut+1:lenCut+Nfe/2,:);
    kspace = rfft*ifftshift(fft(fftshift(imTemp2, 1), [], 1), 1);
end

if ( length(s) == 3 )
    kspaceTemp = fftshift(kspace_in,1);
    imTemp = rifft*ifftshift(ifft(kspaceTemp, [], 1), 1);

    lenCut = Nfe / 4;    
    imTemp2 = imTemp(lenCut+1:lenCut+Nfe/2,:,:);
    kspace = rfft*ifftshift(fft(fftshift(imTemp2, 1), [], 1), 1);
end

if ( length(s) == 4 )
    kspaceTemp = fftshift(kspace_in,1);
    imTemp = rifft*ifftshift(ifft(kspaceTemp, [], 1), 1);

    lenCut = Nfe / 4;    
    imTemp2 = imTemp(lenCut+1:lenCut+Nfe/2,:,:,:);
    kspace = rfft*ifftshift(fft(fftshift(imTemp2, 1), [], 1), 1);
end

if ( length(s) == 5 )
    kspaceTemp = fftshift(kspace_in,1);
    imTemp = rifft*ifftshift(ifft(kspaceTemp, [], 1), 1);

    lenCut = Nfe / 4;    
    imTemp2 = imTemp(lenCut+1:lenCut+Nfe/2,:,:,:,:);
    kspace = rfft*ifftshift(fft(fftshift(imTemp2, 1), [], 1), 1);
end

if ( length(s) == 6 )
    kspaceTemp = fftshift(kspace_in,1);
    imTemp = rifft*ifftshift(ifft(kspaceTemp, [], 1), 1);

    lenCut = Nfe / 4;    
    imTemp2 = imTemp(lenCut+1:lenCut+Nfe/2,:,:,:,:,:);
    kspace = rfft*ifftshift(fft(fftshift(imTemp2, 1), [], 1), 1);
end

if ( length(s) == 7 )
    kspaceTemp = fftshift(kspace_in,1);
    imTemp = rifft*ifftshift(ifft(kspaceTemp, [], 1), 1);

    lenCut = Nfe / 4;    
    imTemp2 = imTemp(lenCut+1:lenCut+Nfe/2,:,:,:,:,:,:);
    kspace = rfft*ifftshift(fft(fftshift(imTemp2, 1), [], 1), 1);
end

if ( length(s) == 8 )
    kspaceTemp = fftshift(kspace_in,1);
    imTemp = rifft*ifftshift(ifft(kspaceTemp, [], 1), 1);

    lenCut = Nfe / 4;    
    imTemp2 = imTemp(lenCut+1:lenCut+Nfe/2,:,:,:,:,:,:,:);
    kspace = rfft*ifftshift(fft(fftshift(imTemp2, 1), [], 1), 1);
end

if ( length(s) == 9 )
    kspaceTemp = fftshift(kspace_in,1);
    imTemp = rifft*ifftshift(ifft(kspaceTemp, [], 1), 1);

    lenCut = Nfe / 4;    
    imTemp2 = imTemp(lenCut+1:lenCut+Nfe/2,:,:,:,:,:,:,:,:);
    kspace = rfft*ifftshift(fft(fftshift(imTemp2, 1), [], 1), 1);
end

if ( length(s) == 10 )
    kspaceTemp = fftshift(kspace_in,1);
    imTemp = rifft*ifftshift(ifft(kspaceTemp, [], 1), 1);

    lenCut = Nfe / 4;    
    imTemp2 = imTemp(lenCut+1:lenCut+Nfe/2,:,:,:,:,:,:,:,:,:);
    kspace = rfft*ifftshift(fft(fftshift(imTemp2, 1), [], 1), 1);
end

if ( length(s) == 11 )
    kspaceTemp = fftshift(kspace_in,1);
    imTemp = rifft*ifftshift(ifft(kspaceTemp, [], 1), 1);

    lenCut = Nfe / 4;    
    imTemp2 = imTemp(lenCut+1:lenCut+Nfe/2,:,:,:,:,:,:,:,:,:,:);
    kspace = rfft*ifftshift(fft(fftshift(imTemp2, 1), [], 1), 1);
end
