function ftkPerfusionDenoisingRun(imageName, filteredImageName, smoothingMode, mode, tuning, neighborHoodSize, tuningFactor, BatFileName, exePath)

if ( isempty(exePath) == 0 )
    command = [exePath 'ftkPerfusionDenoising' ' ' imageName ' ' filteredImageName];
else
    command = ['ftkPerfusionDenoising' ' ' imageName ' ' filteredImageName];
end

%% parameters
if ( isempty(mode) == 0 )
    command = [command ' ' mode];
else
    command = [command ' ' 'WAV'];
end

if ( isempty(smoothingMode) == 0 )
    command = [command ' ' '-smoothing ' smoothingMode];
end

if ( isempty(tuning) == 0 )
    command = [command ' ' '-tuning ' tuning];
else
    command = [command ' ' '-tuning medium'];
end

if ( isempty(neighborHoodSize) == 0 )
    command = [command ' ' '-neighborhood ' num2str(neighborHoodSize)];
end

if ( isempty(tuningFactor) == 0 )
    command = [command ' ' '-tuningFactor ' num2str(tuningFactor)];
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