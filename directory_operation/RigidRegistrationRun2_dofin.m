
function RigidRegistrationRun2_dofin(target, source, dofname, parameterfile, dofinfile)
% rigidly register the image from t2 to tof

if ( isempty(parameterfile) == 0 )
    command = ['rreg' ' ' target ' ' source ' ' '-parameter' ' ' parameterfile ' '  '-dofout' ' ' dofname ' ' '-dofin' ' ' dofinfile];
else
    command = ['rreg' ' ' target ' ' source ' ' '-dofout' ' ' dofname ' ' '-dofin' ' ' dofinfile];
end
disp(command)
[s, w] = dos(command, '-echo');
return