
function TransformationRun(target, source, output, dof)

if (isempty(dir(output))==0)
    disp([' exists : ' output ' ... ' ]);
    return;
end

command = ['transformation' ' ' source ' ' output ' ' '-dofin' ' ' dof ' ' '-target' ' ' target ' '  '-bspline' ];
command
[s, w] = dos(command, '-echo');
return