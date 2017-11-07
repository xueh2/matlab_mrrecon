
function NonRigidRegistrationRun2_2D(target, source, dofnames, parameterfile, numofparameters,...
    dofin, TpValue, controlPoints)
% rigidly register the image from t2 to tof

command = ['hreg2D' ' ' target ' ' source ' ' num2str(numofparameters)];
if ( isempty(parameterfile) == 0 )
    
    command = [command ' '  '-parameter' ];
    
    for i = 1:numofparameters
        command = [command ' ' parameterfile{i}];
    end
end

if ( isempty(dofin) == 0 )
    command = [command ' ' '-dofin' ' ' dofin];
end

command = [command ' '  '-dofout' ];
for i = 1:numofparameters
    command = [command ' ' dofnames{i}];
end
if ( numofparameters == 0 )
    command = [command ' ' dofnames{i}];
end

if ( isempty(TpValue) == 0 )
    command = [command ' ' '-Tp' ' ' num2str(TpValue)];
end

if ( isempty(controlPoints) == 0 )
    command = [command ' ' '-dx' ' ' num2str(controlPoints(1)) ' ' '-dy' ' ' num2str(controlPoints(2)) ' ' '-dz' ' ' num2str(controlPoints(3))];
end

disp(command)

[s, w] = dos(command, '-echo');
return