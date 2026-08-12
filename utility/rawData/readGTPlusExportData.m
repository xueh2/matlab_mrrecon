
function data = readGTPlusExportData(baseFileName)
% data = readGTPlusExportData(baseFileName)

try
    rRes = niftiread([baseFileName '_REAL.hdr']); %NDim = numel(size(rRes)); rRes = permute(rRes, [2 1 3:NDim]);
    iRes = niftiread([baseFileName '_IMAG.hdr']); %NDim = numel(size(iRes)); iRes = permute(iRes, [2 1 3:NDim]);
    data = complex(rRes, iRes);
catch
    data = niftiread([baseFileName '.hdr']);
end
size(data);

