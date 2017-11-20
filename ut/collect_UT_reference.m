function collect_UT_reference(UTCollectRef, resDir, refDir)
% collect the reference

if ( UTCollectRef )
    if ( ~exist(refDir) )
        mkdir(refDir);
    end

    delete(fullfile(refDir, '*.nii'));
    delete(fullfile(refDir, '*.hdr'));
    delete(fullfile(refDir, '*.img'));
    delete(fullfile(refDir, '*.attrib'));
    delete(fullfile(refDir, '*.xml'));
    delete(fullfile(refDir, '*.h5'));
    delete(fullfile(refDir, '*.dcm'));
    delete(fullfile(refDir, '*.Icehead'));
    delete(fullfile(refDir, '*.ima'));
    
    [names, num] = findFILE(resDir, '*.hdr');
    for ii=1:num
        [path, filename, ext] = fileparts(names{ii});
        copyfile(names{ii}, fullfile(refDir, [filename ext]) );
    end

    [names, num] = findFILE(resDir, '*.img');
    for ii=1:num
        [path, filename, ext] = fileparts(names{ii});
        copyfile(names{ii}, fullfile(refDir, [filename ext]) );
    end

    [names, num] = findFILE(resDir, '*.attrib');
    for ii=1:num
        [path, filename, ext] = fileparts(names{ii});
        copyfile(names{ii}, fullfile(refDir, [filename ext]) );
    end
    
    [names, num] = findFILE(resDir, '*.dcm');
    for ii=1:num
        [path, filename, ext] = fileparts(names{ii});
        copyfile(names{ii}, fullfile(refDir, [filename ext]) );
    end
    
    [names, num] = findFILE(resDir, '*.Icehead');
    for ii=1:num
        [path, filename, ext] = fileparts(names{ii});
        copyfile(names{ii}, fullfile(refDir, [filename ext]) );
    end
    
    [names, num] = findFILE(resDir, '*.h5');
    for ii=1:num
        [path, filename, ext] = fileparts(names{ii});
        copyfile(names{ii}, fullfile(refDir, [filename ext]) );
    end
end