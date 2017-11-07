
function ftkDeformationFieldJacobian3DRun(dx, dy, dz, jac, BatFileName, exePath)

if ( isempty(exePath) )    
    command = ['ftkDeformationFieldJacobian3D' '  ' dx ' ' dy ' ' dz ' ' jac];
else
    command = [exePath 'ftkDeformationFieldJacobian3D' '  ' dx ' ' dy ' ' dz ' ' jac ];       
end

disp(command)

if ( isempty(BatFileName) )
    [s, w] = dos(command, '-echo');
else
    fp = fopen(BatFileName, 'a');
    fprintf(fp, '\n');
    fprintf(fp, '%s\n', command);
    fclose(fp);
end
return
