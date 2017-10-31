function ftkMRParameterMappingWithMocoRun(imageName, path, prefix, ... 
    moco, strategy, initial, inv, iters, sigma, neighbor, stepDiv, moreIterInv, algo, ...
    maxIter, FixedDropForInitial, numOfDroppedFrames, framesUsedForIntialFitting, ...
    holeFilling, mapType, algoType, reweighting, costType, T2Prep, T1Inversion, T2Star, ...
    maskType, maskInd, decreasingIntensityCheck, initialMoCoWithRejection, BatFileName, exePath)

if ( isempty(exePath) == 0 )
    command = [exePath 'ftkMRParameterMapping' ' ' '1' ' ' imageName ' ' path ' '];
else
    command = ['ftkMRParameterMapping' ' ' '1' ' ' imageName ' ' path ' '];
end

if ( isempty(prefix) == 0 )
    command = [command ' ' '-prefix' ' ' prefix];
end

%% general moco
if ( isempty(moco) == 0 )
    if( moco )
        command = [command ' ' '-moco' ' ' 'true'];
    else
        command = [command ' ' '-moco' ' ' 'false'];
    end
end

if ( isempty(strategy) == 0 )
    command = [command ' ' '-strategy' ' ' strategy];
end

if ( isempty(initial) == 0 )
    if( initial )
        command = [command ' ' '-initial' ' ' 'true'];
    else
        command = [command ' ' '-initial' ' ' 'false'];
    end
end

if ( isempty(inv) == 0 )
    if( inv )
        command = [command ' ' '-inv' ' ' 'true'];
    else
        command = [command ' ' '-inv' ' ' 'false'];
    end
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

if ( isempty(stepDiv) == 0 )
    command = [command ' ' '-stepDiv' ' ' num2str(stepDiv)];
end

if ( isempty(moreIterInv) == 0 )
    if( moreIterInv )
        command = [command ' ' '-moreIterInv' ' ' 'true'];
    else
        command = [command ' ' '-moreIterInv' ' ' 'false'];
    end
end

if ( isempty(algo) == 0 )
    command = [command ' ' '-algo' ' ' algo];
end

%% model based moco
if ( isempty(maxIter) == 0 )
    command = [command ' ' '-maxIter' ' ' num2str(maxIter)];
end

if ( isempty(FixedDropForInitial) == 0 )
    if( FixedDropForInitial )
        command = [command ' ' '-FixedDropForInitial' ' ' 'true'];
    else
        command = [command ' ' '-FixedDropForInitial' ' ' 'false'];
    end
end

if ( isempty(numOfDroppedFrames) == 0 )    
    command = [command ' ' '-numOfDroppedFrames' ' ' num2str(numel(numOfDroppedFrames)) ' ' num2str(numOfDroppedFrames)];    
end

if ( isempty(framesUsedForIntialFitting) == 0 )    
    command = [command ' ' '-framesUsedForIntialFitting' ' ' num2str(numel(framesUsedForIntialFitting)) ' ' num2str(framesUsedForIntialFitting)];    
end

%% parameter mapping
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

if ( isempty(maskInd) == 0 )
    command = [command ' ' '-maskInd' ' ' num2str(maskInd)];
end

if ( decreasingIntensityCheck )
    command = [command ' ' '-decreasingIntensityCheck true'];
else
    command = [command ' ' '-decreasingIntensityCheck false'];
end

if ( initialMoCoWithRejection )
    command = [command ' ' '-initialMoCoWithRejection true'];
else
    command = [command ' ' '-initialMoCoWithRejection false'];
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