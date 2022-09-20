function copy_data_from_gt_azure_to_azure(dataDir, siteDir, copy_script)
% copy data to azure
% copy_data_from_gt_azure_to_azure('/mnt/Lab-Kellman/RawData', 'BARTS', '/home/xueh2/copy_from_gt_to_barts.sh')

[subdirs, num] = FindSubDirs(fullfile(dataDir, siteDir));

fid = fopen(copy_script, 'w');


% fprintf(fid, '%s\n', 'gt_sas_token="sp=racwdli&st=2022-05-21T14:33:51Z&se=2023-05-21T22:33:51Z&spr=https&sv=2020-08-04&sr=c&sig=TBY6rEBak7mR8iXmHVE61kEfcWJA72W0JOffMCM3aRs%3D"');
% 
% for i=1:num
%     copy_str = ['azcopy copy ' dataDir '/' siteDir '/' subdirs{i} ' https://gadgetronrawdata.blob.core.windows.net/rawdata/' siteDir '?${gt_sas_token} --recursive --block-blob-tier="archive"'];
%     copy_str
%     fprintf(fid, '%s\n', ['echo ' subdirs{i}]);
%     fprintf(fid, '%s\n', copy_str);
% end

fprintf(fid, '%s\n', 'gt_sas_token_src="sp=racwdli&st=2022-05-21T14:33:51Z&se=2023-05-21T22:33:51Z&spr=https&sv=2020-08-04&sr=c&sig=TBY6rEBak7mR8iXmHVE61kEfcWJA72W0JOffMCM3aRs%3D"');
fprintf(fid, '%s\n', 'gt_sas_token_dst="sp=racwdli&st=2022-07-06T14:46:18Z&se=2023-07-06T22:46:18Z&spr=https&sv=2021-06-08&sr=c&sig=Z6PpjssrEXke%2BzZ4oSiSDZ7O9RtZv2b19P1hveTATzE%3D"');

for i=1:num
    copy_str = ['azcopy copy https://gadgetronrawdata.blob.core.windows.net/rawdata/' siteDir '/' subdirs{i} '?${gt_sas_token_src} ' ' https://cmrrawdata.blob.core.windows.net/rawdata/' siteDir '?${gt_sas_token_dst} --recursive --block-blob-tier="archive"'];
    copy_str
    fprintf(fid, '%s\n', ['echo ' subdirs{i}]);
    fprintf(fid, '%s\n', copy_str);
end

fclose(fid);