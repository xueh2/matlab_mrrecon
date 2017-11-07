
function Save3DTo2D_Resize(home, SingleSliceDir, Image3D, perfix, suffix, minW, minH)
% save 3D image to 2D slices

if ( isempty(dir(Image3D)) == 1 )
    return;
end

[data, header] = LoadAnalyze(Image3D, 'Grey');
[pathstr, name, ext, versn] = fileparts(Image3D);

if ( header.xsize<minH & header.ysize<minW )
    dstfile_hres = fullfile(pathstr, [name '_hres' ext]);
    
    ratio = (minH+2)/header.xsize;
    if ( ratio < (minW+2)/header.ysize )
        ratio = (minW+2)/header.ysize;
    end
    
    newXVoxelSize = header.xvoxelsize / ratio;
    newYVoxelSize = header.yvoxelsize / ratio;    
    
    command = ['resample ' Image3D ' ' dstfile_hres ...
        ' -size ' num2str(newXVoxelSize) ...
        ' ' num2str(newYVoxelSize) ...
        ' ' num2str(header.zvoxelsize) ...
        ' -bspline'];
    dos(command, '-echo');
    
    dstfile = fullfile(home, 'ROI', [name ext]);
    movefile(Image3D, dstfile);
    
    srcfile = fullfile(home, [name '.img']);
    dstfile = fullfile(home, 'ROI', [name '.img']);
    movefile(srcfile, dstfile);
    
    [data, header] = LoadAnalyze(dstfile_hres, 'Grey');

end

num = header.zsize;

header2D = header;
header2D.zsize = 1;

for i=1:num
    
    sliceData = data(:, :, i);
    sliceName = fullfile(home, SingleSliceDir, [perfix num2str(i) suffix] );
    SaveAnalyze(sliceData, header2D, sliceName, 'Grey');
end