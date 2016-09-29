
function sampledRangePE = detectSampledRangePE(underSampledKspace)
% detect the sampled data range for the PE direction
% sampledRangePE is [indOfFirstLine indOfLastLine]

sumKspace = underSampledKspace;
for dim=4:9
    sumKspace = sum(sumKspace, dim);
end

loc = detectSampledLines(sumKspace(:,:,:,1));
diff = loc([2:end 1]) - loc;
ind = find(diff==1);
sampledRangePE(1) = loc(ind(1));
sampledRangePE(2) = loc(ind(end))+1;