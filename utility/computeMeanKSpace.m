function meanKspace = computeMeanKSpace(kspace)
% kspace: [COL LIN CHA PHS]
% meanKspace: [COL LIN CHA]

sampledLineLoc = detectSampledLinesIrregular(kspace);

COL = size(kspace, 1);
LIN = size(kspace, 2);
CHA = size(kspace, 3);
PHS = size(kspace, 4);

count = zeros(1, LIN);

for phs=1:PHS
    ind = sampledLineLoc{phs};
    count(ind(:)) = count(ind(:)) + 1;
end

count(find(count==0)) = 1;

countAll = repmat(count, [COL 1 CHA]);

meanKspace = sum(kspace, 4);
meanKspace = meanKspace ./ countAll;
