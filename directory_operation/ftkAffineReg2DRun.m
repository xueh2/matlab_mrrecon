
function ftkAffineReg2DRun(target, source, dofIn, dofOut, parameterfile, transformName, BatFileName)
% rigidly register the image from t2 to tof

command = ['ftkAffineReg2D' ' ' target ' ' source ' ' '-dofout' ' ' dofOut];

if ( isempty(parameterfile) == 0 )
    command = [command ' ' '-parin' ' ' parameterfile];
end

if ( isempty(transformName) == 0 )
    command = [command ' ' '-transform' ' ' transformName];
end

if ( isempty(dofIn) == 0 )
    command = [command ' ' '-dofin' ' ' dofIn];
end

if ( isempty(BatFileName) )
    [s, w] = dos(command, '-echo');
else
    fp = fopen(BatFileName, 'a');
    fprintf(fp, '\n');
    fprintf(fp, '%s\n', command);
    fclose(fp);
end
return