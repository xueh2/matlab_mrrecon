function [flow, meanVelocity, maxVelocity, cardiacoutput] = computeFlow(mag, phs, venc, meanRR, offset, FOV, mask)
% [flow, meanVelocity, maxVelocity] = computeFlow(mag, phs, venc, meanRR, offset, FOV, mask)
% compute flow in ml/second, meanVelocity (cm/s), maxVelocity (cm/s), cardiacoutput (L/min)
% venc: cm/s; meanRR in second

data  = double(mag);
% velocity = double(phs) - offset; velocity = 0.5* venc * velocity/offset;
velocity_rad = pi * (double(phs) - offset)/offset;
velocity = venc * (velocity_rad/pi);

% compute the flow

flow = zeros(size(data, 3), 1); % cm3/s or ml/s
meanVelocity = zeros(size(data, 3), 1); % cm3/s or ml/s
maxVelocity = zeros(size(data, 3), 1); % cm3/s or ml/s

RO = size(data, 1);
E1 = size(data, 1);

pixelArea = ( FOV(1)/RO * FOV(2)/E1 ) / (100); % cm2

for i=1:size(data, 3)
    
    mask2D = mask(:,:,i);
    w = sum(mask2D(:));
    f = mask2D .* velocity(:,:,i);
    
    meanVelocity(i) = sum(f(:))/w;
    maxVelocity(i) = max(f(:));
    
    f = f * pixelArea;
    flow(i) = sum(f(:));       
end

% compute cardiac output
PHS = size(data, 3);
tt = [0 0.5:1:PHS PHS] * (meanRR/PHS);
flowPadded = [flow(1); flow; flow(PHS)];

cardiacoutput = trapz(tt, flowPadded); % ml per meanRR
cardiacoutput = cardiacoutput * (60/meanRR); % ml per minute
cardiacoutput = cardiacoutput / 1e3; % L per second
