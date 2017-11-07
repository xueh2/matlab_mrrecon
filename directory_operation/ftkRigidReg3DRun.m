
function ftkRigidReg3DRun(target, source, dofIn, dofOut, parameterfile, transformName, BatchFile, exePath)
% rigidly register the image

%command = ['ftkRigidReg3D' ' ' target ' ' source ' ' '-dofout' ' ' dofOut];
command = [fullfile(exePath, 'ftkRigidReg3D') ' ' target ' ' source ' ' '-dofout' ' ' dofOut];

if ( isempty(parameterfile) == 0 )
    command = [command ' ' '-parin' ' ' parameterfile];
end

if ( isempty(transformName) == 0 )
    command = [command ' ' '-transform' ' ' transformName];
end

if ( isempty(dofIn) == 0 )
    command = [command ' ' '-dofin' ' ' dofIn];
end

disp(command);

if ( isempty(BatchFile) )
    [s, w] = dos(command, '-echo');
else
    fp = fopen(BatchFile, 'a');
    fprintf(fp, '\n');
    fprintf(fp, '%s\n', command);
    fclose(fp);
end
return