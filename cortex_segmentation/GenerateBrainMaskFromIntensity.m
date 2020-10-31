
function mask = GenerateBrainMaskFromIntensity(filename, maskName)

[data, header] = LoadAnalyze(filename, 'Grey');

mask = zeros(size(data));

mask(find(data>0)) = 1;

SaveAnalyze(uint32(mask), header, maskName, 'Grey');

return