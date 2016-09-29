function CompareTwoGtResults_WholeDirectory(home, folder1, series1, folder2, series2, time_dim)
% CompareTwoGtResults_WholeDirectory(home, folder1, series1, folder2, series2, time_dim)
% CompareTwoGtResults_WholeDirectory('F:\data\CV_Sparse', 'slep_res', 0, 'nl_spirit_res', 1, 3)
% [RO E1 SLC E2 CON PHS REP SET AVE RUN]

[subdirs, num] = FindSubDirs(home);

for n=1:num
       
    subdirs{n}
    a = readGTPlusExportImageSeries(fullfile(home, subdirs{n}, folder1), series1);
    b = readGTPlusExportImageSeries(fullfile(home, subdirs{n}, folder2), series2);
    size(a)
    size(b)
    
    ndim = numel(size(a));
    
    b = b * norm(a(:))/norm(b(:));
    c = cat(3, a, b);
    c = squeeze(c);
    size(c)
    
    if(isempty(time_dim))
        figure; imagescn(c);
    else
        figure; imagescn(c, [], [], [], time_dim);
    end
    
    pause
    closeall    
end