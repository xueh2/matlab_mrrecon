
function ftkRigidReg2DRun(target, source, dofIn, dofOut, parameterfile, transformName, BatchFile)
% rigidly register the image from t2 to tof

%command = ['ftkRigidReg2D' ' ' target ' ' source ' ' '-dofout' ' ' dofOut];
command = ['C:\cc_views\hxue_snapshot_view\SCR_MR_CV\MrFtk\prod\bin\vc9\Release\ftkRigidReg2D' ' ' target ' ' source ' ' '-dofout' ' ' dofOut];

if ( isempty(parameterfile) == 0 )
    command = [command ' ' '-parin' ' ' parameterfile];
end

if ( isempty(transformName) == 0 )
    command = [command ' ' '-transform' ' ' transformName];
end

if ( isempty(dofIn) == 0 )
    command = [command ' ' '-dofin' ' ' dofIn];
end

if ( isempty(BatchFile) )
    [s, w] = dos(command, '-echo');
else
    fp = fopen(BatchFile, 'a');
    fprintf(fp, '\n');
    fprintf(fp, '%s\n', command);
    fclose(fp);
end
return