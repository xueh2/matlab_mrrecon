
function ftkFFDReg2DRun(target, source, numofparameters, dofin, dofnames,...
    parameterfile, controlPoints, transformedFile, dxFile, dyFile, BatFileName)
% rigidly register the image from t2 to tof

command = ['ftkFFDReg2D' '  ' target '  ' source '  ' num2str(numofparameters)];
if ( isempty(parameterfile) == 0 )
    
    command = [command '  '  '-parin' ];
    
    for i = 1:numofparameters
        command = [command '  ' parameterfile{i}];
    end
end

if ( isempty(dofin) == 0 )
    command = [command '  ' '-dofin' '  ' dofin];
end

command = [command '  '  '-dofout' ];
for i = 1:numofparameters
    command = [command '   ' dofnames{i}];
end

if ( isempty(controlPoints) == 0 )
    command = [command '  ' '-ds' '  ' num2str(controlPoints(1)) '  ' num2str(controlPoints(2)) ];
end

if ( (isempty(dxFile)==0) & (isempty(dyFile)==0) )
    command = [command '  ' '-deform' '  ' dxFile '  ' dyFile ];
end

if ( isempty(transformedFile) == 0 )
    command = [command '  ' '-transform' '  ' transformedFile ];
end

disp(command)

if ( isempty(BatFileName) )
    [s, w] = dos(command, '-echo');
else
    fp = fopen(BatFileName, 'a');
    fprintf(fp, '\n');
    fprintf(fp, '%s\n', command);
    fclose(fp);
end
return