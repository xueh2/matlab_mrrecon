function [linearity, LPSF, BWRatio] = quantifyLPSF(LPSFIm, temporalScalingFactorSPIRiT, pts, perturbV)

numOfReg = numel(temporalScalingFactorSPIRiT);
numOfPt = size(pts, 1);

linearity = zeros(numOfPt, numOfReg);

LPSF = cell(numOfReg, numOfPt, 3);

computeLPSF = 1;
if ( size(LPSFIm,4) < 3 )
    computeLPSF = 0;
end

BWRatio = zeros(numOfPt, numOfReg);

isLinear = 0;
if ( size(LPSFIm,5) == 1 )
    isLinear = 1;
end

for tScaling=1:numel(temporalScalingFactorSPIRiT)
      
    if ( isLinear )
        d1 = LPSFIm(:,:,:,2) - LPSFIm(:,:,:,1);
        if ( computeLPSF )
            d2 = LPSFIm(:,:,:,3) - LPSFIm(:,:,:,1);
        end
    else
        d1 = LPSFIm(:,:,:,2,tScaling) - LPSFIm(:,:,:,1,tScaling);
        if ( computeLPSF )
            d2 = LPSFIm(:,:,:,3,tScaling) - LPSFIm(:,:,:,1,tScaling);
        end        
    end
    
    for pt=1:numOfPt 

        currPt = pts(pt, :);

        fil = ones(5, 1)/5;
        w = 3;

        d1Normed = d1 ./ perturbV(pt, 2);
        if ( computeLPSF ) 
            d2Normed = d2 ./ perturbV(pt, 3);
            linearity(pt, tScaling) = abs(d2Normed(pts(pt, 2), pts(pt, 1), pts(pt, 3)))/abs(d1Normed(pts(pt, 2), pts(pt, 1), pts(pt, 3)));
            linearity(pt, tScaling);
        end      
    
        % compute the resolution along  x and y
        upsamplingRatio = 1;
        mainLobeCutoffRatio = 0.64;
        [resolutionFELinear, maxResponseFE, responseFEInterp, ...
            resolutionPELinear, maxResponsePE, responsePEInterp, ...
            resolutionTimeLinear, maxResponseTime, responseTimeInterp] = computeLocalResolutionFromLocalPSF3D(d1, currPt, upsamplingRatio, mainLobeCutoffRatio);

        LPSF{tScaling, pt, 1} = responseFEInterp;
        LPSF{tScaling, pt, 2} = responsePEInterp;
        LPSF{tScaling, pt, 3} = responseTimeInterp;
        
        if ( isLinear )
            responseTimeInterp2 = [zeros(1000, 1); responseTimeInterp; zeros(1000,1)];    
            s = fft(responseTimeInterp2);
            s = fftshift(s);  
            magS = abs(s);
            magS2 = magS;
            equivalentBW = sum(magS2)/numel(s)/max(magS2(:));
            BWRatio(pt, tScaling) = equivalentBW;
        else
            [v, ind] = max(abs(responseTimeInterp));
            responseTimeInterp(1:ind-w) = 0;
            responseTimeInterp(ind+w,end) = 0;

            responseTimeInterp = [zeros(1000, 1); responseTimeInterp; zeros(1000,1)];    
            s = fft(responseTimeInterp);
            s = fftshift(s);

            magS = abs(s);
            magS2 = conv(magS, fil, 'same');

            equivalentBW = sum(magS2)/numel(s)/max(magS2(:));
            BWRatio(pt,tScaling) = equivalentBW;            
        end        
    end
end
