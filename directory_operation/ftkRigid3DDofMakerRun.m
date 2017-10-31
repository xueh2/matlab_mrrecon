
function ftkRigid3DDofMakerRun(dof, tx, ty, tz, rx, ry, rz, BatFileName)
% generate rigid 3D dof file

command = ['ftkRigid3DDofMaker' '  ' dof ]

if ( isempty(tx) == 0 )
    command = [command '  ' '-Tx' '  ' num2str(tx)];
end

if ( isempty(ty) == 0 )
    command = [command '  ' '-Ty' '  ' num2str(ty)];
end

if ( isempty(tz) == 0 )
    command = [command '  ' '-Tz' '  ' num2str(tz)];
end

if ( isempty(rx) == 0 )
    command = [command '  ' '-Rx' '  ' num2str(rx)];
end

if ( isempty(ry) == 0 )
    command = [command '  ' '-Ry' '  ' num2str(ry)];
end

if ( isempty(rz) == 0 )
    command = [command '  ' '-Rz' '  ' num2str(rz)];
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