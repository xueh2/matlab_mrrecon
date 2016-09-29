
function perf_normalized = Perfusion_PrepareSignal_PDNormalization(perf, perf_PD, FA_perf, FA_PD, filtSize)
% perform the signal normalization using PD images
% perf: [RO E1 REP SLC]
% perf_PD: [RO E1 numPD SLC]

%% compute mean PD

RO = size(perf, 1);
E1 = size(perf, 2);

REP = size(perf, 3);
SLC = size(perf, 4);

meanPD = mean(perf_PD, 3);
meanPD = squeeze(meanPD);

perf_normalized = perf;

% for every slice

for s=1:SLC
    
    PD = meanPD(:,:,s);
    
    % smoothing PD with median filter
    % PD_fil = medfilt2(PD, [filtSize filtSize], 'symmetric');
    
    kspaceFilterRO = generateKSpaceFilter('Hanning', 'Medium', RO, [1 RO], RO/2+1);
    kspaceFilterE1 = generateKSpaceFilter('Hanning', 'Medium', E1, [1 E1], E1/2+1);
    
    PD_fil = performRawDataFilter(fft2c(PD), kspaceFilterRO, kspaceFilterE1);
    PD_fil = abs(ifft2c(PD_fil));
    
    scale = mean(PD(:))/mean(PD_fil(:));    
    PD_fil = PD_fil * scale;
    
    ind = find(PD_fil(:)==0);
    PD_fil(ind(:)) = 1;
    
    % compute siganl norm
    perf_norm_slc = perf(:,:,:,s) ./ repmat(PD_fil, [1 1 REP]);   
    perf_norm_slc = perf_norm_slc * sin(FA_PD*pi/180)/sin(FA_perf*pi/180);
    
    perf_normalized(:,:,:,s) = perf_norm_slc;
end
