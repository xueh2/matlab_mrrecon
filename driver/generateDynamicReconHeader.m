
function header = generateDynamicReconHeader(f, reductionFactor, blockCoord, Nfe, Npe, numOfCoil, minFEUsed, maxFEUsed, sampledLineLoc)
% generate the header structure for frame f in the dynamic recon setting

header = struct('asymm_echo', 0, 'center_kspace_fe', 0, 'skip_points', 0, 'subsampling_factor', 1, ...
    'blocks', 1, 'central_kspace', 0, 'coil_type', 2, 'fft', 2, 'percent', 80, 'Nfe', 0, 'Npe', 0, 'num_coils', 0, ...
    'sampling_location', [], 'consider_coils', 0, 'center_locs', 0);

header.subsampling_factor = reductionFactor;
header.blocks = length(blockCoord); 
header.fft = 2;
header.Nfe = Nfe;
header.Npe = Npe;
header.num_coils = numOfCoil; 
header.minFEUsed = ceil(minFEUsed*Nfe);
header.maxFEUsed = floor(maxFEUsed*Nfe);

header.sampling_location = sampledLineLoc(:,f);
header.grappaACSLines = findGrappaAcsLines(header.sampling_location, reductionFactor, Npe);
