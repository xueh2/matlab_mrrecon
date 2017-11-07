
function CutImageRun(currentDir, tofDir, tofCut, background)

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

[xy, yz, zx]=getMIP3(data, header);


%background = 30;
areaPercentage = 0.5;
[thresholdxy, segmentedMIPxy] = OptimalThreshold(xy, areaPercentage, background);
[thresholdyz, segmentedMIPyz] = OptimalThreshold(yz, areaPercentage, background);
[thresholdzx, segmentedMIPzx] = OptimalThreshold(zx, areaPercentage, background);

t = min([thresholdxy, thresholdyz, thresholdzx]);
t = t * 0.7;

background = t * 0.7;

[seedx, seedy, seedz, numofseeds, largestComponent] = zbsfindseeds_nobackground(data, header, t, background);

%figure; plotSeeds(seedx, seedy, seedz, numofseeds, 0.2, [], 'b.');

filename = fullfile(currentDir, tofCut, 'seeds.mat');
save(filename, 'seedx', 'seedy', 'seedz', 'largestComponent');

[ROI, headerROI, leftup, rightdown] = CutImageArea2(data, header, largestComponent, t);

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