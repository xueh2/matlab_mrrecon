function DicomFixImageOrientation(fname, dst_fname, flipX, flipY)
% DicomFixImageOrientation(fname, dst_fname, flipX, flipY)

if(~flipX & ~flipY) return; end

info = dicominfo(fname);
[x, map] = dicomread(fname);

if(flipX)
    x = flipdim(x, 2);
end

if(flipY)
    x = flipdim(x, 1);
end

if(isempty(map))
    dicomwrite(x, dst_fname, info);
else
    dicomwrite(x, map, dst_fname, info);
end