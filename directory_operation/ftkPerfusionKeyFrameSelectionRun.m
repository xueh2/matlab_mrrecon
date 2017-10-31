function ftkPerfusionKeyFrameSelectionRun(imageName, keyFrame, numOfPDMaps, uniformKeyFrame, maxKeyFrameRatio, BatFileName, exePath)

if ( isempty(exePath) == 0 )
    command = [exePath 'ftkPerfusionKeyFrameSelection' ' ' '1' ' ' imageName ' '];
else
    command = ['ftkPerfusionKeyFrameSelection' ' ' '1' ' ' imageName ' '];
end

if ( isempty(keyFrame) == 0 )
    command = [command ' ' '-keyFrame' ' ' num2str(keyFrame)];
end

if ( isempty(numOfPDMaps) == 0 )
    command = [command ' ' '-numOfPD' ' ' num2str(numOfPDMaps)];
end

if ( isempty(maxKeyFrameRatio) == 0 )
    command = [command ' ' '-maxKeyFrameRatio' ' ' num2str(maxKeyFrameRatio)];
end

if ( isempty(uniformKeyFrame) == 0 )
    if ( uniformKeyFrame > 0 )
        command = [command ' ' '-uniformKeyFrame true'];
    else
        command = [command ' ' '-uniformKeyFrame false'];        
    end
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