function CutImageRun_useT2(currentDir, tofDir, tofCut, warpedMask)

% cut current image and put the results into the tofCut

% load tof image
path_abo = fullfile(currentDir, tofDir, '*.hdr' );
indir = dir(path_abo);
num_indir = length(indir);

if ( num_indir == 0 )
    disp('hdr file is not available...');
    return;
end

filename = fullfile(currentDir, tofDir, indir(1).name);
[data, headerdata] = LoadAnalyze(filename, 'Grey');

% load maskmask
path_abo = fullfile(currentDir, warpedMask, '*.hdr' );
indir = dir(path_abo);
num_indir = length(indir);

if ( num_indir == 0 )
    disp('hdr file is not available...');
    return;
end

filename = fullfile(currentDir, warpedMask, indir(1).name);
[mask, headermask] = LoadAnalyze(filename, 'Grey');

mask = logical(mask);

% mask the brain
tofMask = false(size(data));
tofMask(find(data>0)) = 1;

clear data

newMask = Mask_interSect_DITK(tofMask, headerdata, mask, headermask);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_abo = fullfile(currentDir, tofDir, '*.hdr' );
indir = dir(path_abo);
num_indir = length(indir);

if ( num_indir == 0 )
    disp('hdr file is not available...');
    return;
end

filename = fullfile(currentDir, tofDir, indir(1).name);
[data, headerdata] = LoadAnalyze(filename, 'Grey');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SaveAnalyze(newMask, headerdata, 'newMask.hdr', 'Grey');

[ROI, headerROI, leftup, rightdown] = CutImageArea2(data, headerdata, newMask, 0);

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