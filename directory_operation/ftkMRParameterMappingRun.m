function ftkMRParameterMappingRun(imageName, path, prefix, holeFilling, mapType, algoType, reweighting, costType, T2Prep, T1Inversion, T2Star, maskType, BatFileName, exePath)

if ( isempty(exePath) == 0 )
    command = [exePath 'ftkMRParameterMapping' ' ' '1' ' ' imageName ' ' path ' '];
else
    command = ['ftkMRParameterMapping' ' ' '1' ' ' imageName ' ' path ' '];
end

if ( isempty(prefix) == 0 )
    command = [command ' ' '-prefix' ' ' prefix];
end

if ( isempty(holeFilling) == 0 )
    if( holeFilling )
        command = [command ' ' '-holeFilling' ' ' 'true'];
    else
        command = [command ' ' '-holeFilling' ' ' 'false'];
    end
end

if ( isempty(mapType) == 0 )
    command = [command ' ' '-mapType' ' ' mapType];
end

if ( isempty(algoType) == 0 )
    command = [command ' ' '-algoType' ' ' algoType];
end

if ( isempty(reweighting) == 0 )
    if( reweighting )
        command = [command ' ' '-reweighting' ' ' 'true'];
    else
        command = [command ' ' '-reweighting' ' ' 'false'];
    end
end

if ( isempty(costType) == 0 )
    command = [command ' ' '-costType' ' ' costType];
end

if ( isempty(T2Prep) == 0 )
    command = [command ' ' '-parameter T2Prep' ' ' num2str(numel(T2Prep)) ' ' num2str(T2Prep)];
end

if ( isempty(T1Inversion) == 0 )
    command = [command ' ' '-parameter T1Inversion' ' ' num2str(numel(T1Inversion)) ' ' num2str(T1Inversion')];
end

if ( isempty(T2Star) == 0 )
    command = [command ' ' '-parameter T2Star' ' ' num2str(numel(T2Star)) ' ' num2str(T2Star)];
end

if ( isempty(maskType) == 0 )
    command = [command ' ' '-maskType' ' ' maskType];
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