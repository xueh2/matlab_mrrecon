
function CutImageRun_manual(currentDir, tofDir, tofCut, leftup, rightdown)

% cut current image and put the results into the tofCut

path_abo = fullfile(currentDir, tofDir, '*.hdr' );
indir = dir(path_abo);
num_indir = length(indir);

if ( num_indir == 0 )
    disp('hdr file is not available...');
    return;
end

% load image
filename = fullfile(currentDir, tofDir, indir(1).name);
[data, header] = LoadAnalyze(filename, 'Grey');


[ROI, headerROI] = getROI(data, header, leftup, rightdown);

[ROIxy, ROIyz, ROIzx]=getMIP3(ROI, headerROI);
filename = fullfile(currentDir, tofCut, 'ROIxy.tiff');
imwrite(uint8(ROIxy), filename, 'tiff');
filename = fullfile(currentDir, tofCut, 'ROIyz.tiff');
imwrite(uint8(ROIyz), filename, 'tiff');
filename = fullfile(currentDir, tofCut, 'ROIzx.tiff');
imwrite(uint8(ROIzx), filename, 'tiff');

% save
place = find(currentDir(length(currentDir):-1:1) == filesep );
prefix = currentDir(length(currentDir)-place(1)+2:length(currentDir));

filename = fullfile(currentDir, tofCut, [prefix '_tofCut.mat']);
save(filename, 'ROI', 'headerROI', 'leftup', 'rightdown');

filename = fullfile(currentDir, tofCut, [prefix '_tofCut.hdr']);
SaveAnalyze(ROI, headerROI, filename, 'Grey');

filename = fullfile(currentDir, tofCut, [prefix '_ROI.txt']);
fid = fopen(filename, 'w');
fprintf(fid, '%f  %f  %f \n', leftup(1), leftup(2), leftup(3));
fprintf(fid, '%f  %f  %f \n', rightdown(1), rightdown(2), rightdown(3));
fclose(fid);

disp('CutImageRun finished...');

return