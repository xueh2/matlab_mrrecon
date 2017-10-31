function ftkMRPDECoreMOCO2DRun(imageName, MOCOImage, keyFrame, nonRigid_strategy, nonRigid_initial, iters, sigma, neighbor, numOfPDMaps, inverseFlag, volumePreserving, BatFileName, exePath)
if ( isempty(exePath) == 0 )
    command = [exePath 'ftkMrPDECoreMOCO2D' ' ' '1' ' ' imageName ' '];
else
    command = ['ftkMrPDECoreMOCO2D' ' ' '1' ' ' imageName ' '];
end

if ( isempty(MOCOImage) == 0 )
    command = [command ' ' '-reg' ' ' MOCOImage];
end

if ( isempty(nonRigid_strategy) == 0 )
    command = [command ' ' '-strategy' ' ' nonRigid_strategy];
end

if ( isempty(nonRigid_initial) == 0 )
    if( nonRigid_initial )
        command = [command ' ' '-initial' ' ' 'true'];
    else
        command = [command ' ' '-initial' ' ' 'false'];
    end
end

if ( isempty(inverseFlag) == 0 )
    if( inverseFlag )
        command = [command ' ' '-inv' ' ' 'true'];
    else
        command = [command ' ' '-inv' ' ' 'false'];
    end
end

if ( isempty(keyFrame) == 0 )
    command = [command ' ' '-keyFrame' ' ' num2str(keyFrame)];
end

if ( isempty(numOfPDMaps) == 0 )
    command = [command ' ' '-numOfPre' ' ' num2str(numOfPDMaps)];
end

if ( isempty(iters) == 0 )
    command = [command ' ' '-iters' ' ' num2str(numel(iters)) ' ' num2str(iters)];
end

if ( isempty(sigma) == 0 )
    command = [command ' ' '-sigma' ' ' num2str(sigma)];
end

if ( isempty(neighbor) == 0 )
    command = [command ' ' '-neighbor' ' ' num2str(neighbor)];
end

if ( isempty(volumePreserving) == 0 )
    if( volumePreserving )
        command = [command ' ' '-volumePreserving' ' ' 'true'];
    else
        command = [command ' ' '-volumePreserving' ' ' 'false'];
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