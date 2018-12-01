
function header = CreateFtkHeaderInfo(data, voxelsize)

[ysize, xsize, zsize, tsize, nsize, msize] = size(data);

if ( ~isempty(voxelsize) )
    header = struct('sizeX', xsize, 'sizeY', ysize, 'sizeZ', zsize, 'sizeT', tsize, 'sizeN', nsize, 'sizeM', msize, ... 
        'spacingX', voxelsize(1), 'spacingY', voxelsize(2), 'spacingZ', voxelsize(3), ... 
        'spacingT', 1.0, 'spacingN', 1.0, 'spacingM', 1.0, 'positionPatient', [0 0 0], 'orientationPatient', eye(3,3));
else
    header = struct('sizeX', xsize, 'sizeY', ysize, 'sizeZ', zsize, 'sizeT', tsize, 'sizeN', nsize, 'sizeM', msize, ... 
        'spacingX', 1.0, 'spacingY', 1.0, 'spacingZ', 1.0, ... 
        'spacingT', 1.0, 'spacingN', 1.0, 'spacingM', 1.0, 'positionPatient', [0 0 0], 'orientationPatient', eye(3,3));
end

header.xsize = xsize;
header.ysize = ysize;
header.zsize = zsize;
header.xvoxelsize = voxelsize(1);
header.yvoxelsize = voxelsize(2);
header.zvoxelsize = voxelsize(3);