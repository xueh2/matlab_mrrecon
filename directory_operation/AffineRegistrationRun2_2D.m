
function AffineRegistrationRun2_2D(target, source, dofname, parameterfile, dofin)

command = ['areg2D' ' ' target ' ' source ' '  '-dofout' ' ' dofname];
if ( isempty(parameterfile) == 0 )
    command = [command ' '  '-parameter' ' ' parameterfile];
end

if ( isempty(dofin) == 0 )
    command = [command ' ' '-dofin' ' ' dofin];
end

disp(command)

[s, w] = dos(command, '-echo');
return