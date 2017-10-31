function ftkPerfusionParameterMapRun(imageName, pathForMaps, method, prefix, sigmas, filterNeighbor, thresGradient, BatFileName, exePath)

if ( isempty(exePath) == 1 )
    command = ['ftkPerfusionParameterMap' ' ' ' 1 ' imageName ' ' pathForMaps ' ' '-method ' method];
else
    command = [exePath 'ftkPerfusionParameterMap' ' ' ' 1 ' imageName ' ' pathForMaps ' ' '-method ' method];
end

if ( isempty(prefix) == 0 )
    command = [command ' ' '-prefix' ' ' prefix];
end

if ( isempty(sigmas) == 0 )
    command = [command ' ' '-parameter Sigmas' ' ' num2str(length(sigmas)) ' ' num2str(sigmas)];
end

if ( isempty(filterNeighbor) == 0 )
    command = [command ' ' '-parameter Filtering' ' ' num2str(filterNeighbor)];
end

if ( isempty(thresGradient) == 0 )
    command = [command ' ' '-parameter ThresGradient' ' ' num2str(thresGradient)];
end


%% batch file

disp('----------------------------------------------')
disp(command)
disp('----------------------------------------------')

if ( isempty(BatFileName) )
    [s, w] = dos(command, '-echo');
else
    fp = fopen(BatFileName, 'a');
    fprintf(fp, '\n\n');
    fprintf(fp, '%s\n', command);
    fprintf(fp, '\n\n');
    fclose(fp);
end
return