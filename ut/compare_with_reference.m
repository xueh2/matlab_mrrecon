function compare_with_reference(testInfo, resDir, refDir, logfile)
% compare the results and output UT log

fid = fopen(logfile, 'a+');

[names, num] = findFILE(resDir, '*.hdr');
[names2, num2] = findFILE(refDir, '*.hdr');

use_dicom = 0;
if(num==0)
    [names, num] = findFILE(resDir, '*.dcm');
    [names2, num2] = findFILE(refDir, '*.dcm');    
    use_dicom = 1;
end

fprintf(fid, '\n\n%s\n', '---------------------------------------------------------------------------------------');
fprintf(fid, '%s\n', date);
fprintf(fid, '%s\n', '-----------------------------------------------------------------');
fprintf(fid, '%s -- %s\n', testInfo{1}, testInfo{2});
fprintf(fid, '%s -- %s -- %s\n', testInfo{3}, testInfo{4}, testInfo{5});
fprintf(fid, '%s\n', '-----------------------------------------------------------------');
fprintf(fid, 'This run generates %d images \n', num);
fprintf(fid, 'Ref run generates %d images \n', num2);
fprintf(fid, '%s\n', '-----------------------------------------------------------------');

allRight = 1;
for ii=1:num
    
    [path, filename, ext] = fileparts(names{ii});
    
    try
        if(~use_dicom)
            res = analyze75read(names{ii}); 
        else
            res = dicomread(names{ii}); 
        end
        res = double(res);
    catch
        fprintf(fid, '--> Image cannot be read for %d : %s \n', ii, names{ii});
        allRight = 0;
        continue;
    end
    
    try
        if(~use_dicom)
            ref = analyze75read( fullfile(refDir, [filename ext]) ); 
        else
            ref = dicomread( fullfile(refDir, [filename ext]) ); 
        end
        ref = double(ref);
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
    
    diff = (res-ref);   
    d = norm(diff(:))/(norm(ref(:))+eps);
    
    if ( d > 0.01 )
        fprintf(fid, '--> Image difference (%f) is big for %d : %s \n', d, ii, names{ii});
        allRight = 0;
        warning(['--> Image difference ' num2str(d) ' is big for ' num2str(ii) ' : ' names{ii}]);
    else
        % fprintf(fid, '--> Image difference (%f) is good for %d : %s \n', d, ii, names{ii});
    end
end

if ( allRight )
    fprintf(fid, 'Test Succeeds \n');
end
fprintf(fid, '%s\n', '-----------------------------------------------------------------');

fclose(fid);
