
function TransformationRun_BatFile(target, source, output, dof, BatFileName)

if (isempty(dir(output))==0)
    disp([' exists : ' output ' ... ' ]);
%     return;
end

command = ['transformation' ' ' source ' ' output ' ' '-dofin' ' ' dof ' ' '-target' ' ' target ' '  '-bspline' ];
command

if ( isempty(dir(BatFileName)) )
    [s, w] = dos(command, '-echo');
else
    fp = fopen(BatFileName, 'a');
    fprintf(fp, '\n');
    fprintf(fp, '%s\n', command);
    fclose(fp);
end

return