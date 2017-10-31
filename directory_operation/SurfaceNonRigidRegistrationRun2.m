
function SurfaceNonRigidRegistrationRun2(target, source, dofnames, numofparameters,...
    dofin, locator, iterations, spacing, epsilon, symmetric)
% non-rigidly register the two surfaces
% <-locator>          Locator: 0 = cell locator, 1 = point locator, 2 = kd-tree locator (default = 1)
% <-dofin name>       Name of input file
% <-dofout name>      Name of output file
% <-epsilon>          Value for espilon (default=0.01)
% <-symmetric>        Use symmetric distance (default OFF)
% <-ignoreedges>      Ignores edges in ICP (default OFF)
% <-iterations>       Number of 3D registration iterations (default 100)
% <-ds spacing>       Control point spacing

if ( isempty(dir(dofnames{numofparameters}))==0 )
    disp([dofnames{numofparameters} ' already exists ...']);
    return;
end

command = ['shreg_multipleDOFs' ' ' target ' ' source ' ' num2str(numofparameters)];

if ( isempty(dofin) == 0 )
    command = [command ' ' '-dofin' ' ' dofin];
end

command = [command ' '  '-dofout' ];

if ( numofparameters > 1 )
    for i = 1:numofparameters
        command = [command ' ' dofnames{i}];
    end
else
    command = [command ' ' dofnames];
end
% if ( numofparameters == 0 )
%     command = [command ' ' dofnames{i}];
% end

if ( isempty(locator) == 0 )
    command = [command ' ' '-locator' ' ' num2str(locator)];
end

if ( isempty(iterations) == 0 )
    command = [command ' ' '-iterations' ' ' num2str(iterations)];
end

if ( isempty(spacing) == 0 )
    command = [command ' ' '-ds' ' ' num2str(spacing)];
end

if ( isempty(epsilon) == 0 )
    command = [command ' ' '-epsilon' ' ' num2str(epsilon)];
end

if ( isempty(symmetric) == 0 )
    command = [command ' ' '-symmetric'];
end

disp(command)

[s, w] = dos(command, '-echo');
return


