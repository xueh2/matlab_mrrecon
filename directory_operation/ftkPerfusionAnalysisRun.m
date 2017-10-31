function ftkPerfusionAnalysisRun(imageName, ... 
    heartDetection, croppedheart, roiFile, keyFrame, ...
    rigid, rigid_reg, rigid_strategy, rigid_initial, rigid_parin, ...
    affine, affine_reg, affine_strategy, affine_initial, affine_parin, ...
    nonRigid, nonRigid_reg, nonRigid_strategy, nonRigid_initial, iters, sigma, neighbor, ...
    temporalSmoothing, temporal_result, smoothingStrategy, smoothingTuning,  ...
    spatialSmoothing, spatial_result, neighborhood,  ...                                
    scaleSpace, scaleSpacePrefixes, sampleInterval, scaleSigmas, thresSigma, thresGradient, ...
    fullyQuantitative, fullyPrefixes, ...
    deltaT, numOfPDMaps, BatFileName)

command = ['ftkPerfusionAnalysis' ' ' '1' ' ' imageName ' '];

%% heart detection
if ( isempty(heartDetection) == 0 )
    if( heartDetection )
        command = [command ' ' '-heart-detection' ' ' 'true'];
    else
        command = [command ' ' '-heart-detection' ' ' 'false'];
    end
end

if ( isempty(croppedheart) == 0 )
    command = [command ' ' '-croppedheart' ' ' croppedheart];
end

if ( isempty(roiFile) == 0 )
    command = [command ' ' '-roi' ' ' roiFile];
end

if ( isempty(keyFrame) == 0 )
    command = [command ' ' '-keyframe' ' ' num2str(keyFrame)];
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

%% temporal smoothing
if ( isempty(temporalSmoothing) == 0 )
    if( temporalSmoothing )
        command = [command ' ' '-temporal-smoothing' ' ' 'true'];
    else
        command = [command ' ' '-temporal-smoothing' ' ' 'false'];
    end
end

if ( isempty(temporal_result) == 0 )
    command = [command ' ' '-temporal' ' ' temporal_result];
end

if ( isempty(smoothingStrategy) == 0 )
    command = [command ' ' '-smoothing-strategy' ' ' smoothingStrategy];
end

if ( isempty(smoothingTuning) == 0 )
    command = [command ' ' '-smoothing-tuning' ' ' smoothingTuning];
end

%% spatial smoothing
if ( isempty(spatialSmoothing) == 0 )
    if( spatialSmoothing )
        command = [command ' ' '-spatial-smoothing' ' ' 'true'];
    else
        command = [command ' ' '-spatial-smoothing' ' ' 'false'];
    end
end

if ( isempty(spatial_result) == 0 )
    command = [command ' ' '-spatial' ' ' spatial_result];
end

if ( isempty(neighborhood) == 0 )
    command = [command ' ' '-neighborhood' ' ' num2str(neighborhood)];
end

%% semi map generation
if ( isempty(scaleSpace) == 0 )
    if( scaleSpace )
        command = [command ' ' '-scale-space' ' ' 'true'];
    else
        command = [command ' ' '-scale-space' ' ' 'false'];
    end
end

if ( isempty(scaleSpacePrefixes) == 0 )
    command = [command ' ' '-scale-space-prefixes' ' ' scaleSpacePrefixes];
end

if ( isempty(sampleInterval) == 0 )
    command = [command ' ' '-sampleInterval' ' ' num2str(sampleInterval)];
end

if ( isempty(scaleSigmas) == 0 )
    command = [command ' ' '-scaleSigmas' ' ' num2str(numel(scaleSigmas)) ' ' num2str(scaleSigmas)];
end

if ( isempty(thresSigma) == 0 )
    command = [command ' ' '-thresSigma' ' ' num2str(thresSigma)];
end

if ( isempty(thresGradient) == 0 )
    command = [command ' ' '-thresGradient' ' ' num2str(thresGradient)];
end

%% fully perfusion map generation
if ( isempty(fullyQuantitative) == 0 )
    if( fullyQuantitative )
        command = [command ' ' '-fully-quantitative' ' ' 'true'];
    else
        command = [command ' ' '-fully-quantitative' ' ' 'false'];
    end
end

if ( isempty(scaleSpacePrefixes) == 0 )
    command = [command ' ' '-fully-prefixes' ' ' fullyPrefixes];
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