
function sampledLineLoc = detectSampledLines(kspace)
% find the sampled lines
% kspace [fe * pe * coil]

Npe = size(kspace, 2);
sampledLineLoc = [];
for pe = 1:Npe
    data = sum(kspace(:, pe, :));
    if ( sum(abs(data)) > 0 )
        sampledLineLoc = [sampledLineLoc; pe];
    end
end