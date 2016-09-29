
function imPF = PerformPartialFourierHandling2D(im, option, verbose)
% im: complex image [COL LIN CHA N]
% option.rangeUsedFE: the sampled range, FE direction
% option.rangeUsedLIN: the sampled range, FE direction
% option.PFMethod:  'zero-filling', 'homodyne', 'FengHuang', 'SPIRiTPFCorr'
% 'zero-filling': zero-filling the region
% 'homodyne': iterative homodyne filter
% option.homodyne.iters: number of iterations
% option.homodyne.thres: threshold to stop the iteration
% 'FengHuang': if CHA==1, single channel FengHuang method; if CHA>1,
% multi-channel kernel is estimated
% 'SPIRiTPFCorr': only available for multi-channel cases

imPF = im;

if ( isempty(option.rangeUsedFE) & isempty(option.rangeUsedLIN) )    
    return;
end

COL = size(im, 1);
LIN = size(im, 2);
CHA = size(im, 3);
N = size(im, 4);

startCOL = 1;
endCOL = COL;

startLIN = 1;
endLIN = LIN;

if ( ~isempty(option.rangeUsedFE) )
    startCOL = option.rangeUsedFE(1);
    endCOL = option.rangeUsedFE(2);
end

if ( ~isempty(option.rangeUsedLIN) )
    startLIN = option.rangeUsedLIN(1);
    endLIN = option.rangeUsedLIN(2);
end

if ( strcmp(option.PFMethod, 'zero-filling') )
    kspace = fft2c(im);
    kspace(startCOL:endCOL, :) = 0;
    kspace(:, startLIN:endLIN) = 0;
    imPF = ifft2c(kspace);
end

if ( strcmp(option.PFMethod, 'homodyne') )

    % plotComplexImage(im);
    
    imPF = im;
    
    filterType = 'Hanning';
    filterStrength = 'Medium';
    
    if ( isempty(option.rangeUsedFE) )
        filterCOL = generateKSpaceFilter(filterType, filterStrength, COL, [1 COL], floor(COL/2));    
    else
        range = findSymmetricSampledRegion(option.rangeUsedFE, floor(COL/2)+1);
        filterCOL = generateKSpaceFilter(filterType, filterStrength, COL, range, floor(COL/2));
    end
    
    if ( isempty(option.rangeUsedLIN) )
        filterLIN = generateKSpaceFilter(filterType, filterStrength, LIN, [1 LIN], floor(LIN/2));
    else
        filterLIN = generateKSpaceFilter(filterType, filterStrength, LIN, option.rangeUsedLIN, floor(LIN/2));
    end
    
    niter = 1;
    kspace = fft2c(imPF);
    kspaceOri = kspace;
    imPrev = imPF;
    
    for niter=1:option.homodyne.iters
        kspaceFiltered = performRawDataFilter(kspace, filterCOL, filterLIN);
        imFiltered = ifft2c(kspaceFiltered);
        backgroundPhs = imFiltered;
        backgroundPhs = backgroundPhs ./ (abs(backgroundPhs)+eps);
        imPF = imPF .* conj(backgroundPhs);
        kspace = fft2c(imPF);
        
        thresCurr = norm(imPrev(:)-imPF(:))/norm(imPF(:));
        if ( verbose )
            disp(['Homodyne iter : ' num2str(niter) ' - thres : ' num2str(thresCurr) ' ... ']);
        end
        
        if ( thresCurr < option.homodyne.thres )
            break;
        end
        
        imPrev = imPF;
        
        % plotComplexImage(imPF);
    end
    
    rangeFE = [];
    if ( isempty(option.rangeUsedFE) )
        rangeFE = [1 COL];
    else
        rangeFE = option.rangeUsedFE;
    end
    
    rangeLIN = [];
    if ( isempty(option.rangeUsedLIN) )
        rangeLIN = [1 LIN];
    else
        rangeLIN = option.rangeUsedLIN;
    end
    
    kspace(rangeFE(1):rangeFE(2), rangeLIN(1):rangeLIN(2), :, :) = kspaceOri(rangeFE(1):rangeFE(2), rangeLIN(1):rangeLIN(2), :, :);
    imPF = ifft2c(kspace);
end

if ( strcmp(option.PFMethod, 'pocs') )

    % plotComplexImage(im);
    
    imPF = im;
    
    filterType = 'Hanning';
    filterStrength = 'Medium';
    
    rangeFE = [];
    if ( isempty(option.rangeUsedFE) )
        filterCOL = generateKSpaceFilter(filterType, filterStrength, COL, [1 COL], floor(COL/2));    
        rangeFE = [1 COL];
    else
        range = findSymmetricSampledRegion(option.rangeUsedFE, floor(COL/2)+1);
        filterCOL = generateKSpaceFilter(filterType, filterStrength, COL, range, floor(COL/2));
        rangeFE = option.rangeUsedFE;
    end
    
    rangeLIN = [];
    if ( isempty(option.rangeUsedLIN) )
        filterLIN = generateKSpaceFilter(filterType, filterStrength, LIN, [1 LIN], floor(LIN/2));
        rangeLIN = [1 LIN];
    else
        filterLIN = generateKSpaceFilter(filterType, filterStrength, LIN, option.rangeUsedLIN, floor(LIN/2));
        rangeLIN = option.rangeUsedLIN;
    end
    
    kspace = fft2c(imPF);
    kspaceOri = kspace;
    kspaceFiltered = performRawDataFilter(kspace, filterCOL, filterLIN);
    imFiltered = ifft2c(kspaceFiltered);
    phi = imFiltered ./ (abs(imFiltered) + eps);
    imPrev = imPF;
    
    for niter=1:option.pocs.iters
                
        imPocs = abs(imPF).*phi;
        kspacePocs = fft2c(imPocs);
        % plotKSpaceSamplingPattern(kspacePocs);
        
        kspacePocs(rangeFE(1):rangeFE(2), rangeLIN(1):rangeLIN(2), :, :) = kspaceOri(rangeFE(1):rangeFE(2), rangeLIN(1):rangeLIN(2), :, :); 
        
        imPF = ifft2c(kspacePocs);
        
        thresCurr = norm(imPrev(:)-imPF(:))/norm(imPF(:));
        if ( verbose )
            disp(['Pocs iter : ' num2str(niter) ' - thres : ' num2str(thresCurr) ' ... ']);
        end
        
        if ( thresCurr < option.pocs.thres )
            break;
        end
        
        imPrev = imPF;
        
        % plotComplexImage(imPF);
    end
end

if ( strcmp(option.PFMethod, 'FengHuang') )
    kspace = fft2c(im);
    
    rangeFE = [];
    if ( isempty(option.rangeUsedFE) )
        rangeFE = [1 COL];
    else
        rangeFE = option.rangeUsedFE;
    end
    
    rangeLIN = [];
    if ( isempty(option.rangeUsedLIN) )
        rangeLIN = [1 LIN];
    else
        rangeLIN = option.rangeUsedLIN;
    end

    kspaceConj = complexConjugateParitalFourier2D(kspace, rangeFE, rangeLIN);
    % plotKSpaceSamplingPattern(kspace);
    % plotKSpaceSamplingPattern(kspaceConj);

    % region for kernel estimation
    
    rangeSymmetricFE = findSymmetricSampledRegion(rangeFE, floor(COL/2)+1);
    rangeSymmetricLIN = findSymmetricSampledRegion(rangeLIN, floor(LIN/2)+1);
    
    src = kspaceConj(rangeSymmetricFE(1):rangeSymmetricFE(2), rangeSymmetricLIN(1):rangeSymmetricLIN(2), :, :);
    dst = kspace(rangeSymmetricFE(1):rangeSymmetricFE(2), rangeSymmetricLIN(1):rangeSymmetricLIN(2), :, :);
    % plotKSpaceSamplingPattern(src);
    % plotKSpaceSamplingPattern(dst);
    
    src = mean(src, 4);
    dst = mean(dst, 4);
    
    [kernel, rawkernel] = PFConvKernelCalibration(src, dst, option.FengHuang.kSize, option.FengHuang.thresReg);    
    kspacePF = PFConvRecon2D(kspaceConj, kspace, rangeFE, rangeLIN, option.FengHuang.kSize, kernel, rawkernel);
    % plotKSpaceSamplingPattern(kspacePF);
    
    imPF = ifft2c(kspacePF);
end

if ( strcmp(option.PFMethod, 'SPIRiTPFCorr') )
    kspace = fft2c(im);
    
    rangeFE = [];
    if ( isempty(option.rangeUsedFE) )
        rangeFE = [1 COL];
    else
        rangeFE = option.rangeUsedFE;
    end
    
    rangeLIN = [];
    if ( isempty(option.rangeUsedLIN) )
        rangeLIN = [1 LIN];
    else
        rangeLIN = option.rangeUsedLIN;
    end

    kspaceConj = complexConjugateParitalFourier2D(kspace, rangeFE, rangeLIN);
    % plotKSpaceSamplingPattern(kspace);
    % plotKSpaceSamplingPattern(kspaceConj);

    % region for kernel estimation
    
    rangeSymmetricFE = findSymmetricSampledRegion(rangeFE, floor(COL/2)+1);
    rangeSymmetricLIN = findSymmetricSampledRegion(rangeLIN, floor(LIN/2)+1);
    
    src = kspace(rangeSymmetricFE(1):rangeSymmetricFE(2), rangeSymmetricLIN(1):rangeSymmetricLIN(2), :, :);
    dst = kspace(rangeSymmetricFE(1):rangeSymmetricFE(2), rangeSymmetricLIN(1):rangeSymmetricLIN(2), :, :);
    % plotKSpaceSamplingPattern(src);
    % plotKSpaceSamplingPattern(dst);
    
    src = mean(src, 4);
    dst = mean(dst, 4);
    
    headerSrc = CreateFtkHeaderInfo(src, [1 1 1 1]);
    headerDst = CreateFtkHeaderInfo(dst, [1 1 1 1]);

    [kernelS2D, kernelD2S, kernelImS2D, conjKernelImS2D, kernelImD2S, conjKernelImD2S, GP_Im, identityKernelIm, GP_I_Im, conjGP_I_Im, identityKernel] = Matlab_PerformSPIRiTKernelCalibrationSrcDst2D(single(src), headerSrc, single(dst), headerDst, ... 
        option.SPIRiTPFCorr.kSize, option.SPIRiTPFCorr.OutPattern, [LIN COL], option.SPIRiTPFCorr.thresReg);
    
    header = CreateFtkHeaderInfo(kspace, [1 1 1 1]);
    headerKernel = CreateFtkHeaderInfo(GP_I_Im, [1 1 1 1]);            

    kspaceDst = zeros(size(kspace));
    kspaceDst(rangeFE(1):rangeFE(2), rangeLIN(1):rangeLIN(2), :, :) = kspace(rangeFE(1):rangeFE(2), rangeLIN(1):rangeLIN(2), :, :);
    
    kspaceInitial = kspaceConj;
    kspaceInitial(rangeFE(1):rangeFE(2), rangeLIN(1):rangeLIN(2), :, :) = kspace(rangeFE(1):rangeFE(2), rangeLIN(1):rangeLIN(2), :, :);
    
    kspacePF = Matlab_PerformSPIRiTLinearRecon2D(single(kspaceDst), header, single(GP_I_Im), headerKernel, single(kspaceInitial), option.SPIRiTPFCorr.performCenteredFFT, option.SPIRiTPFCorr.maxIter, option.SPIRiTPFCorr.stopThres, option.SPIRiTPFCorr.printInfo);    
            
    % plotKSpaceSamplingPattern(kspacePF);
    
    imPF = ifft2c(kspacePF);
end
