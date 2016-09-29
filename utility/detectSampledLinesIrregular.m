
function sampledLineLoc = detectSampledLinesIrregular(kspace)
% find the sampled lines
% kspace [fe * pe * coil * frame]

Nframe = size(kspace, 4);
sampledLineLoc = cell(Nframe, 1);
for f=1:Nframe
    loc = detectSampledLines(kspace(:,:,:,f));
    sampledLineLoc{f} = loc;
end
