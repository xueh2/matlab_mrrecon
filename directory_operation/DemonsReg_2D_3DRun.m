
function DemonsReg_2D_3DRun(target, index, parameterfile, regVolume, dx_file, dy_file)
% rigidly register the image from t2 to tof

command = ['demons_insights_2D_3D'  '  '  target  '  '  index  '  '  regVolume '  '  '-parin'  ' '  parameterfile '  '  ...
        '-DeformationField'  '  '  dx_file '  ' dy_file];

disp(command)
[s, w] = dos(command, '-echo');

command = ['transformation_field_2D_3D' ' '  target  ' ' regVolume ' ' '-Physical' ' '  '-DeformationField'  '  '  dx_file '  ' dy_file ' '  '-bspline'];
disp(command)
[s, w] = dos(command, '-echo');

return