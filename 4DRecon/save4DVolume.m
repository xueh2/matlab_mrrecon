function save4DVolume(data4D, header, dstFolder, prefix, sizeRatio, centre, width, delay)

if ( ~exist('prefix') )
    prefix = 'Data4D';
end

if ( ~exist('sizeRatio') )
    sizeRatio = 2;
end

if ( ~exist('centre') )
    centre = 300;
end

if ( ~exist('prefix') )
    width = 500;
end

if ( ~exist('delay') )
    delay = 1/30;
end

NSlc = size(data4D, 3);
NPhs = size(data4D, 4);

for s=1:NSlc    
    [im3D, headerSlc] = getSLC(data4D, header, s);    
    filename = fullfile(dstFolder, [prefix '_SLC' num2str(s)]);
    Matlab_SaveAnalyze(single(im3D), headerSlc, [filename '.hdr']);
    hdr2gifWithWindowSetting([filename '.hdr'], [filename '.gif'], sizeRatio, centre, width, delay);
end

for b=1:NPhs   
    [im3D, headerPhs] = getPHS(data4D, header, b);
    filename = fullfile(dstFolder, [prefix '_PHS' num2str(b)]);
    Matlab_SaveAnalyze(single(im3D), headerPhs, [filename '.hdr']);
    hdr2gifWithWindowSetting([filename '.hdr'], [filename '.gif'], sizeRatio, centre, width, delay);
end
