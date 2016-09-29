function copyAnalzyeFile(srcDir, dstDir, srcFileName, dstFileName)

if ( isFileExist(fullfile(srcDir, [srcFileName '.hdr'])) )
    copyfile(fullfile(srcDir, [srcFileName '.hdr']), fullfile(dstDir, [dstFileName '.hdr']));
    copyfile(fullfile(srcDir, [srcFileName '.img']), fullfile(dstDir, [dstFileName '.img']));
end