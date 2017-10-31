
function TransformationDeformFieldsRun_BatFile2(target, source, output,...
    dof, dxFileName, dyFileName, dzFileName, BatFileName, Interpolation)

% if (isempty(dir(output))==0)
%     disp([' exists : ' output ' ... ' ]);
% %     return;
% end

command = ['transformation_DeformFields' ' ' source ' ' output ' ' '-dofin' ' ' dof...
    ' ' Interpolation ];

if (isempty(target)==0)
    command = [command ' ' '-target ' target];
end

if (isempty(dxFileName)==0)
    command = [command ' ' '-DeformX ' dxFileName];
end

if (isempty(dyFileName)==0)
    command = [command ' ' '-DeformY ' dyFileName];
end

if (isempty(dzFileName)==0)
    command = [command ' ' '-DeformZ ' dzFileName];
end

command

if ( isempty(BatFileName) )
    BatFileName = '';
end

if ( isempty(dir(BatFileName)) )
    [s, w] = dos(command, '-echo');
else
    fp = fopen(BatFileName, 'a');
    fprintf(fp, '\n');
    fprintf(fp, '%s\n', command);
    fclose(fp);
end

return