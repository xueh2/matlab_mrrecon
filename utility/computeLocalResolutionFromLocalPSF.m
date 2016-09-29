function [resolutionFE, maxResponseFE, responseFEInterp, resolutionPE, maxResponsePE, responsePEInterp] = computeLocalResolutionFromLocalPSF(LPSF, upsamplingRatio, mainLobeCutoffRatio)
% compute the resolution from the local point spread function

Nfe = size(LPSF, 1);
Npe = size(LPSF, 2);

LPSFMag = abs(LPSF);
[I, J] = find(LPSFMag==max(LPSFMag(:)));

responseFE = squeeze(LPSF(:, J));
responseFE = reshape(responseFE, [Nfe 1]);
responsePE = squeeze(LPSF(I, :));
responsePE = reshape(responsePE, [Npe 1]);

responseFEInterp = Zero_Padding_Resize_NoFiltering(responseFE,upsamplingRatio*Nfe, 1);
responsePEInterp = Zero_Padding_Resize_NoFiltering(responsePE,upsamplingRatio*Npe, 1);

% responseFEInterp = interp(abs(responseFE),upsamplingRatio);
% responsePEInterp = interp(abs(responsePE),upsamplingRatio);

% plot(abs(responseFEInterp));
% plot(abs(responsePEInterp));

% ratio = upsamplingRatio;
% LPSFInterp = Zero_Padding_Resize_NoFiltering(LPSF,ratio*Nfe,ratio*Npe);
% LPSFInterpMag = abs(LPSFInterp);
% [I, J] = find(LPSFInterpMag==max(LPSFInterpMag(:)));

% FE
response = responseFEInterp;
[maxResponseFE, peakInd] = max(response);

locR = peakInd;
for loc=peakInd:size(response,1)
    if ( response(loc) <= maxResponseFE*mainLobeCutoffRatio )
        locR = loc;
        break;
    end
end

locL = peakInd;
for loc=peakInd:-1:1
    if ( response(loc) <= maxResponseFE*mainLobeCutoffRatio )
        locL = loc;
        break;
    end
end

resolutionFE = (locR-locL+1)/upsamplingRatio;

% PE
response = responsePEInterp;
[maxResponsePE, peakInd] = max(response);

locR = peakInd;
for loc=peakInd:size(response,1)
    if ( response(loc) <= maxResponsePE*mainLobeCutoffRatio )
        locR = loc;
        break;
    end
end

locL = peakInd;
for loc=peakInd:-1:1
    if ( response(loc) <= maxResponsePE*mainLobeCutoffRatio )
        locL = loc;
        break;
    end
end

resolutionPE = (locR-locL+1)/upsamplingRatio;
