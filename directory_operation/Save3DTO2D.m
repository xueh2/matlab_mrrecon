
function Save3DTo2D(home, SingleSliceDir, Image3D, perfix, suffix)
% save 3D image to 2D slices

if ( isempty(dir(Image3D)) == 1 )
    return;
end

[data, header] = Matlab_LoadAnalyze(Image3D);
num = header.zsize;

header2D = header;
header2D.zsize = 1;
header2D.sizeZ = 1;

for i=1:num
    
    sliceData = data(:, :, i);
    sliceName = fullfile(home, SingleSliceDir, [perfix num2str(i) suffix] );
    Matlab_SaveAnalyze(int16(sliceData), header2D, sliceName);
end