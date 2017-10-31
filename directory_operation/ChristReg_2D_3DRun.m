
function ChristReg_2D_3DRun(target, index, parameterfile, regVolume, dx_file, dy_file)
% rigidly register the image from t2 to tof

command = ['ChristReg_2D_3D_Smoothed'  '  '  target  '  '  index  '  '  regVolume '  '  '-parin'  ' '  parameterfile '  '  '-bspline'  '  '...
        '-DeformationField'  '  '  dx_file '  ' dy_file];

disp(command)
[s, w] = dos(command, '-echo');
return