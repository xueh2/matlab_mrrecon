function data = removeOutlier(data, percentageOutlier)
% remove large outliear

[nn, v] = hist(data(:), 1024);

nn = nn(end:-1:1);
nn = cumsum(nn);
nn = nn ./ numel(data);
p = find(nn<=percentageOutlier);

if ( ~isempty(p) )
    outlierInd = numel(nn) - p(end) + 1;
    data(find(data>=v(outlierInd))) = v(outlierInd)-1;
end