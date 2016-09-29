
function sampledLineLoc = detectSampledLinesDynamic(kspace)
% find the sampled lines
% kspace [fe * pe * coil * frame]

Nframe = size(kspace, 4);
sampledLineLoc = [];
for f=1:Nframe
    loc = detectSampledLines(kspace(:,:,:,f));
    if ( isempty(loc) )
        continue;
    end
    
    if ( ~isempty(sampledLineLoc) )
        if ( numel(loc) >= size(sampledLineLoc, 1) )

            R = loc(2)-loc(1);
            supposedLines = loc(1):R:loc(end);

            if ( numel(loc) == numel(supposedLines) )
                sampledLineLoc = [sampledLineLoc loc];
            else
                sampledLineLoc = [sampledLineLoc supposedLines'];
            end
        else
            maxLoc = max(sampledLineLoc(:));
            redutionFactor = sampledLineLoc(2,1)-sampledLineLoc(1,1);
            loc = loc(1):redutionFactor:maxLoc;
            loc = loc';
            if ( numel(loc) ~= size(sampledLineLoc, 1) )
                loc = [loc; loc(end)+redutionFactor];
                sampledLineLoc = [sampledLineLoc loc];
            else
                sampledLineLoc = [sampledLineLoc loc];
            end
        end
    else
        sampledLineLoc = [sampledLineLoc loc];
    end
end
