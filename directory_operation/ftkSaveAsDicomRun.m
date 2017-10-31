function ftkSaveAsDicomRun(imageName, resultDir, templateDir, seriesNumber, interactiveFlag, BatFileName, exePath)

if ( isempty(exePath) == 0 )
    command = [exePath 'ftkSaveAsDicom' ' ' imageName ' ' resultDir];
else
    command = ['ftkSaveAsDicom' ' ' imageName ' ' resultDir];
end

%% parameters
if ( isempty(templateDir) == 0 )
    command = [command ' ' '-database ' templateDir];
end

if ( isempty(seriesNumber) == 0 )
    command = [command ' ' '-series ' num2str(seriesNumber)];
else
    command = [command ' ' '-series 2009'];
end

if ( interactiveFlag )
    command = [command ' ' '-interactive true'];
else
    command = [command ' ' '-interactive false'];
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