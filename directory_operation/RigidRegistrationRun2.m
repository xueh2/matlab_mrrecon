
function RigidRegistrationRun2(target, source, dofname, parameterfile)
% rigidly register the image from t2 to tof

if ( isempty(parameterfile) == 0 )
    command = ['rreg' ' ' target ' ' source ' ' '-parameter' ' ' parameterfile ' '  '-dofout' ' ' dofname];
else
    command = ['rreg' ' ' target ' ' source ' ' '-dofout' ' ' dofname];
end
disp(command)
[s, w] = dos(command, '-echo');
return