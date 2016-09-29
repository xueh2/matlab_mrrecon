function compare_with_reference_ICE(testInfo, resDir, refDir, logfile)
% compare the results and output UT log

fid = fopen(logfile, 'a+');

[names, num] = findFILE(resDir, '*.dcm');
[names2, num2] = findFILE(refDir, '*.dcm');

fprintf(fid, '\n\n%s\n', '---------------------------------------------------------------------------------------');
fprintf(fid, '%s\n', date);
fprintf(fid, '%s\n', '-----------------------------------------------------------------');
fprintf(fid, '%s -- %s\n', testInfo{1}, testInfo{2});
fprintf(fid, '%s -- %s \n', testInfo{3}, testInfo{4});
fprintf(fid, '%s\n', '-----------------------------------------------------------------');
fprintf(fid, 'This run generates %d images \n', num);
fprintf(fid, 'Ref run generates %d images \n', num2);
fprintf(fid, '%s\n', '-----------------------------------------------------------------');

allRight = 1;
for ii=1:num
    
    [path, filename, ext] = fileparts(names{ii});
    
    try
        res = dicomread(names{ii}); res = double(res);
    catch
        fprintf(fid, '--> Image cannot be read for %d : %s \n', ii, names{ii});
        allRight = 0;
        continue;
    end
    
    try
        ref = dicomread( fullfile(refDir, [filename ext]) ); ref = double(ref);
    catch
        fprintf(fid, '--> Image cannot be read for %d : %s \n', ii, names2{ii});
        allRight = 0;
        continue;
    end
    
    s = size(res);
    s2 = size(ref);
    
    if ( s(1)~=s2(1) | s(2)~=s2(2) | s(2)~=s2(2) )
        fprintf(fid, '--> Image size does not match for %d : %s \n', ii, names{ii});
        allRight = 0;
        continue;
    end
    
    diff = res-ref;
    
    dRef = norm(ref(:));
    d = norm(diff(:));
    
    if ( d/dRef > 0.75 )
        fprintf(fid, '--> Image difference (%f) is big for %d : %s \n', d/dRef, ii, names{ii});
        allRight = 0;
        warning(['--> Image difference ' num2str(d/dRef) ' is big for ' num2str(ii) ' : ' names{ii}]);
    end
end

if ( allRight )
    fprintf(fid, 'Test Succeeds \n');
end
fprintf(fid, '%s\n', '-----------------------------------------------------------------');

fclose(fid);
