function copy_data_to_azure(dataDir, siteDir, copy_script)
% copy data to azure

[subdirs, num] = FindSubDirs(fullfile(dataDir, siteDir));

fid = fopen(copy_script, 'w');


fprintf(fid, '%s\n', 'gt_sas_token="sp=racwdli&st=2022-05-21T14:33:51Z&se=2023-05-21T22:33:51Z&spr=https&sv=2020-08-04&sr=c&sig=TBY6rEBak7mR8iXmHVE61kEfcWJA72W0JOffMCM3aRs%3D"');

for i=1:num
    copy_str = ['azcopy copy ' dataDir '/' siteDir '/' subdirs{i} ' https://gadgetronrawdata.blob.core.windows.net/rawdata/' siteDir '?${gt_sas_token} --recursive --block-blob-tier="archive"'];
    copy_str
    fprintf(fid, '%s\n', ['echo ' subdirs{i}]);
    fprintf(fid, '%s\n', copy_str);
end

fclose(fid);