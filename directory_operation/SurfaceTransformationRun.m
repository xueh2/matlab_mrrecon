
function SurfaceTransformationRun(source, output, dof, invert)

if ( isempty(dir(output))==0 )
    disp([output ' already exists ...']);
    return;
end

command = ['stransformation' ' ' source ' ' output ' ' '-dofin' ' ' dof ];

if ( exist('invert') == 1 )
    command = [command ' ' '-invert'];
end

command
[s, w] = dos(command, '-echo');
return