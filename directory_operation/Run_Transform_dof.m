
function Run_Transform_dof(SourceFile, doffile, warpedFile, warpedBrainMask)

% transform the image
command = ['transformation' ' ' SourceFile ' ' 'test.hdr' ' ' '-dofin' ' ' doffile ' ' '-bspline'];

disp(command)

[s, w] = dos(command, '-echo');

[data, header] = LoadAnalyze('test.hdr', 'Grey');
se = strel('ball', 3, 3, 0);
pp = imclose(data, se);
SaveAnalyze(uint32(pp), header, 'test.hdr', 'Grey');
        

command = ['subtract' ' ' warpedBrainMask ' ' 'test.hdr' ' ' warpedFile ' ' '-no_norm'];
disp(command)
[s, w] = dos(command, '-echo');

[data, header] = LoadAnalyze(warpedFile, 'Grey');        
[T, header] = LoadAnalyze(warpedBrainMask, 'Grey');        

tt = T - data;
SaveAnalyze(uint32(tt), header, warpedFile, 'Grey');

delete('test.hdr');
delete('test.img');
return