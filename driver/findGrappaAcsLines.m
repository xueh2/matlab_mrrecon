
function grappaACSLines = findGrappaAcsLines(sampling_location, reductionFactor, Npe)
% find the ACS lines

senLines = setdiff(sampling_location, 1:reductionFactor:Npe);
if ( isempty(senLines) )
    grappaACSLines = [];
else
    grappaACSLines = senLines(1)-1:senLines(end)+1;
end
