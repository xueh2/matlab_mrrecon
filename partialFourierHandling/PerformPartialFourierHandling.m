
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

if ( empty(option.rangeUsedFE) & empty(option.rangeUsedLIN) )
    imPF = im;
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

if ( ~empty(option.rangeUsedFE) )
    startCOL = option.rangeUsedFE(1);
    endCOL = option.rangeUsedFE(2);
end

if ( ~empty(option.rangeUsedLIN) )
    startLIN = option.rangeUsedLIN(1);
    endLIN = option.rangeUsedLIN(2);
end

kspace = fft2c(im);

if ( strcmp(option.PFMethod, 'zero-filling') )
    kspace(startCOL:endCOL, :) = 0;
    kspace(:, startLIN:endLIN) = 0;
    imPF = ifft2c(kspace);
end

if ( strcmp(option.PFMethod, 'homodyne') )

    imPF = im;
    
    filterType = 'Hanning';
    filterStrength = 'Medium';
    
    filterCOL = generateKSpaceFilter(filterType, filterStrength, COL, option.rangeUsedFE, floor(COL/2));
    filterLIN = generateKSpaceFilter(filterType, filterStrength, LIN, option.rangeUsedLIN, floor(LIN/2));
    
    niter = 1;
    kspace = fft2c(imPF);
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
    end
end
