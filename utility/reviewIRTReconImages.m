function vol = reviewIRTReconImages(resFolder, startInd, endInd, prefix)
% review the IRF results

if ( nargin < 4 )
    prefix = 'IRT_REP';
end

filename = fullfile(resFolder, [prefix '_' num2str(startInd) '_mag.hdr']);    
[data, header] = Matlab_LoadAnalyze(filename);

vol = zeros(size(data, 1), size(data, 2), endInd-startInd+1);

for ii=startInd:endInd    
    filename = fullfile(resFolder, [prefix '_' num2str(ii) '_mag.hdr']);    
    [data, header] = Matlab_LoadAnalyze(filename);
    vol(:,:,ii) = data;
end

figure; imagescn(vol, [], [], [], 3);