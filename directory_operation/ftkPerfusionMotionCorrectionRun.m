function ftkPerfusionMotionCorrectionRun(imageName, mocStrategy, ... 
    keyFrame, regionofheart, roiFile, croppedheart, ...
    rigid, rigid_reg, rigid_strategy, rigid_initial, rigid_parin, ...
    affine, affine_reg, affine_strategy, affine_initial, affine_parin, ...
    nonRigid, nonRigid_reg, nonRigid_strategy, nonRigid_initial, iters, sigma, neighbor, ...
    deltaT, numOfPDMaps, config, BatFileName, exePath)

if ( isempty(exePath) == 0 )
    command = [exePath 'ftkPerfusionMotionCorrection' ' ' '1' ' ' imageName ' '];
else
    command = ['ftkPerfusionMotionCorrection' ' ' '1' ' ' imageName ' '];
end

%% heart detection
if ( isempty(mocStrategy) == 0 )
    command = [command ' ' '-strategy' ' ' mocStrategy];
end

if ( isempty(keyFrame) == 0 )
    command = [command ' ' '-keyframe' ' ' num2str(keyFrame)];
end

if ( isempty(regionofheart) == 0 )
    command = [command ' ' '-regionofheart' ' ' num2str(regionofheart)];
end

if ( isempty(roiFile) == 0 )
    command = [command ' ' '-roi' ' ' roiFile];
end

if ( isempty(croppedheart) == 0 )
    command = [command ' ' '-croppedheart' ' ' croppedheart];
end

%% rigid
if ( isempty(rigid) == 0 )
    if( rigid )
        command = [command ' ' '-rigid' ' ' 'true'];
    else
        command = [command ' ' '-rigid' ' ' 'false'];
    end
end

if ( isempty(rigid_reg) == 0 )
    command = [command ' ' '-rigid_reg' ' ' rigid_reg];
end

if ( isempty(rigid_strategy) == 0 )
    command = [command ' ' '-rigid_strategy' ' ' rigid_strategy];
end

if ( isempty(rigid_initial) == 0 )
    if( rigid_initial )
        command = [command ' ' '-rigid_initial' ' ' 'true'];
    else
        command = [command ' ' '-rigid_initial' ' ' 'false'];
    end
end

if ( isempty(rigid_parin) == 0 )
    command = [command ' ' '-rigid_parin' ' ' rigid_parin];
end

%% affine
if ( isempty(affine) == 0 )
    if( affine )
        command = [command ' ' '-affine' ' ' 'true'];
    else
        command = [command ' ' '-affine' ' ' 'false'];
    end
end

if ( isempty(affine_reg) == 0 )
    command = [command ' ' '-affine_reg' ' ' affine_reg];
end

if ( isempty(affine_strategy) == 0 )
    command = [command ' ' '-affine_strategy' ' ' affine_strategy];
end

if ( isempty(affine_initial) == 0 )
    if( affine_initial )
        command = [command ' ' '-affine_initial' ' ' 'true'];
    else
        command = [command ' ' '-affine_initial' ' ' 'false'];
    end
end

if ( isempty(affine_parin) == 0 )
    command = [command ' ' '-affine_parin' ' ' affine_parin];
end

%% non-rigid

if ( isempty(nonRigid) == 0 )
    if( nonRigid )
        command = [command ' ' '-non-rigid' ' ' 'true'];
    else
        command = [command ' ' '-non-rigid' ' ' 'false'];
    end
end

if ( isempty(nonRigid_reg) == 0 )
    command = [command ' ' '-non-rigid_reg' ' ' nonRigid_reg];
end

if ( isempty(nonRigid_strategy) == 0 )
    command = [command ' ' '-non-rigid_strategy' ' ' nonRigid_strategy];
end

if ( isempty(nonRigid_initial) == 0 )
    if( nonRigid_initial )
        command = [command ' ' '-non-rigid_initial' ' ' 'true'];
    else
        command = [command ' ' '-non-rigid_initial' ' ' 'false'];
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

%% config
if ( isempty(config) == 0 )
    command = [command ' ' '-config' ' ' config];
end

%% data properties
if ( isempty(deltaT) == 0 )
    command = [command ' ' '-deltaT' ' ' num2str(deltaT)];
end

if ( isempty(numOfPDMaps) == 0 )
    command = [command ' ' '-numOfPDMaps' ' ' num2str(numOfPDMaps)];
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