
function ind = findHighProb(prob, probablity_thres, minN)
% find the index of items with high probabilities

ind = find(prob>=probablity_thres);

while ( length(ind)<minN )
    probablity_thres = probablity_thres - 0.0005;
    if ( probablity_thres <= 0.02 )
        ind = find(prob>=probablity_thres);
        break;
    end
    ind = find(prob>=probablity_thres);
end

if ( length(ind) == 0 )
    warning(' length(ind) == 0 ...');
end

return
