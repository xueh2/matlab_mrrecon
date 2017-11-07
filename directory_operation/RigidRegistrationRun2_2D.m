
function RigidRegistrationRun2_2D(target, source, dofname, parameterfile)
% rigidly register the image from t2 to tof

if ( isempty(parameterfile) == 0 )
    command = ['rreg2D' ' ' target ' ' source ' ' '-parameter' ' ' parameterfile ' '  '-dofout' ' ' dofname];
else
    command = ['rreg2D' ' ' target ' ' source ' ' '-dofout' ' ' dofname];
end
disp(command)
[s, w] = dos(command, '-echo');
return