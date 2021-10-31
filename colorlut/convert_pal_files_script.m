% script to convert Numaris4 color LUT pal files to NumarisX pal files

% pal file directories
idir = 'Z:\Share\MATLAB\MrRecon\colorlut\ColorLUT_VE11E';
odir = 'Z:\Share\MATLAB\MrRecon\colorlut\ColorLUT_XA20A';

idir = 'Z:\Share\MATLAB\MrRecon\colorlut\ColorLUT_VE11E';
odir = 'Z:\Share\MATLAB\MrRecon\colorlut\ColorLUT_XA20A';

d = dir([idir,filesep,'*pal']);
for i = 1:length(d)
    disp(['reading: ',idir,filesep,d(i).name])
    [header,rgb]=read_color_lut_pal_file([idir,filesep,d(i).name]);
    disp(['writing: ',odir,filesep,d(i).name])
    write_color_lut_NX([odir,filesep,d(i).name],rgb)
end

return