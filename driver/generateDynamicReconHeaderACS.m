
function headerAcs = generateDynamicReconHeaderACS(header, reconStrategy, linesUsedForGRAPPA, reductionFactor, Npe)
% generate header for acs data

headerAcs = header;
if ( strcmp(reconStrategy, lower('GRAPPA')) | strcmp(reconStrategy, lower('SENSE')) )

    if ( strcmp(linesUsedForGRAPPA, lower('frame')) )
        firstAcsLine = header.grappaACSLines(1);
        lastAcsLine = header.grappaACSLines(end);
        senLines = setdiff(header.grappaACSLines - firstAcsLine:reductionFactor:lastAcsLine);
        senLines = senLines';
        minLine = min(header.sampling_location);
        while ( senLines(1) < minLine )
            senLines = circshift(senLines, -1);
        end
        headerAcs.sampling_location = senLines;                
    elseif ( strcmp(linesUsedForGRAPPA, lower('All')) )
        firstAcsLine = header.grappaACSLines(1);
        lastAcsLine = header.grappaACSLines(end);
        senLines = setdiff(header.grappaACSLines, firstAcsLine:reductionFactor:lastAcsLine);
        senLines = senLines';
        minLine = min(header.sampling_location);
        while ( senLines(1) < minLine )
            senLines = circshift(senLines, -1);
        end
        headerAcs.sampling_location = senLines;
        for r=1:reductionFactor-2
            headerAcs.sampling_location = [headerAcs.sampling_location; senLines+r]; 
        end
        headerAcs.sampling_location = mod(headerAcs.sampling_location, Npe);
        headerAcs.sampling_location(find(headerAcs.sampling_location==0)) = Npe;

        sampling_locationAcs = [];
        for l=1:numel(headerAcs.sampling_location)
            if ( ~isempty(find(header.grappaACSLines==headerAcs.sampling_location(l))) )
                sampling_locationAcs = [sampling_locationAcs; headerAcs.sampling_location(l)];
            end
        end
        headerAcs.sampling_location = sampling_locationAcs;
    end
else
    if ( strcmp(linesUsedForGRAPPA, lower('frame')) )
        senLines = setdiff(1:Npe, headerAcs.sampling_location);
        senLines = senLines';
        minLine = min(header.sampling_location);
        while ( senLines(1) < minLine )
            senLines = circshift(senLines, -1);
        end
        headerAcs.sampling_location = senLines;
    elseif ( strcmp(linesUsedForGRAPPA, lower('All')) )
        senLines = setdiff(1:Npe, headerAcs.sampling_location);
        senLines = senLines';
        minLine = min(header.sampling_location);
        while ( senLines(1) < minLine )
            senLines = circshift(senLines, -1);
        end
        headerAcs.sampling_location = senLines;
        for r=1:reductionFactor-2
            headerAcs.sampling_location = [headerAcs.sampling_location; senLines+r]; 
        end
        headerAcs.sampling_location = mod(headerAcs.sampling_location, Npe);
        headerAcs.sampling_location(find(headerAcs.sampling_location==0)) = Npe;
    end
end
    