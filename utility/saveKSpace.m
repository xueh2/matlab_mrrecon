function saveKSpace(kspace, dstDir, prefix, pixelSpacing)
% save kspace as analyze format
% kspace : [COL LIN CHA]

if ( ~exist('pixelSpacing') )
    pixelSpacing = [1 1 1];
end

header = CreateFtkHeaderInfo(kspace, pixelSpacing);

filename = fullfile(dstDir, [prefix '_real.hdr']);
Matlab_SaveAnalyze(single(real(kspace)), header, filename);

filename = fullfile(dstDir, [prefix '_imag.hdr']);
Matlab_SaveAnalyze(single(imag(kspace)), header, filename);

filename = fullfile(dstDir, [prefix '_Mag.hdr']);
Matlab_SaveAnalyze(single(abs(kspace)), header, filename);

filename = fullfile(dstDir, [prefix '_Phase.hdr']);
Matlab_SaveAnalyze(single(angle(kspace)), header, filename);