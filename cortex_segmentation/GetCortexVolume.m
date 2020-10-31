
function cortex_volume = GetCortexVolume(Internal, External, brainmask, header)
% get cortex volume
% the minThickness constraint is used
% cortex_volume: mL, cm^3

index = find( (Internal>=0) & (External<=0) & (brainmask>0) );
num = length(index);
cortex_volume = num*header.xvoxelsize*header.yvoxelsize*header.zvoxelsize* 10^-3;

return