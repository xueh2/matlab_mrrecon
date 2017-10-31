
function Run_Transform_dof_Mask(MaskFile, doffile, ImageFile, WarpedBrainMask, MaskedImage)

% transform the image
command = ['transformation' ' ' MaskFile ' ' WarpedBrainMask ' ' '-dofin' ' ' doffile ' ' '-bspline'];

disp(command)

[s, w] = dos(command, '-echo');

[data, header] = LoadAnalyze(WarpedBrainMask, 'Grey');
se = strel('ball', 3, 3, 0);
pp = imclose(data, se);
SaveAnalyze(uint32(pp), header, WarpedBrainMask, 'Grey');
        
[data, header] = LoadAnalyze(ImageFile, 'Grey');
data(find(pp==0)) = 0;
SaveAnalyze(data, header, MaskedImage, 'Grey');

return