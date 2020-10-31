

% create vtk surfaces from the binary image

inputfile = 'J:\more_neonates_images_LS\templates\complex_brain\seg\deepgray.hdr';
outputfile = 'J:\more_neonates_images_LS\templates\complex_brain\seg\deepgray.vtk';

[data, header] = LoadAnalyze(inputfile, 'Grey');

data(find(data>0)) = 1;
threshold = 0.5;
shrink = [1 1 1];
smooth = [0.1 20];

mcubes2VTKFile(single(data), header, outputfile, threshold, shrink, smooth);
