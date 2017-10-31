
function Transformation_N_Run(target, source, output, dofs, N, BatFileName)

if (isempty(dir(output))==0)
    disp([' exists : ' output ' ... ' ]);
%     return;
end

command = ['transformationN'  ' '  source  ' '  output  ' '  num2str(N)  ' '  '-target'  ' '  target  ' '  '-bspline'  ' '  '-dofin' ];

for i=1:N
    
   command = [command ' ' dofs{i} ' ' ];
    
end

command

if ( isempty(BatFileName) )
    [s, w] = dos(command, '-echo');
else
    fid = fopen(BatFileName, 'a');
    fprintf(fid, '\n%s\n', command);
    fclose(fid);
end

return