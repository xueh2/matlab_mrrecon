function CopyGadgetronInstallersAll(dstDir, chroot_string, gt_string)
% copy gadgetron installer to the dst folder
% CopyGadgetronInstallersAll('\\137.187.134.135\share\Installer', chroot_string, gt_string)

% CopyGadgetronInstallers('D:\installer\VB17', fullfile(dstDir, 'VB17'), chroot_string);
% CopyGadgetronInstallers('D:\installer\VB17_Seq', fullfile(dstDir, 'VB17'), gt_string);

% CopyGadgetronInstallers('D:\installer\VD13', fullfile(dstDir, 'VD13'), chroot_string);
% CopyGadgetronInstallers('D:\installer\VD13_Seq', fullfile(dstDir, 'VD13'), gt_string);

%CopyGadgetronInstallers('D:\installer\VE11', fullfile(dstDir, 'VE11'), chroot_string);
% CopyGadgetronInstallers('D:\installer\VE11_Seq', fullfile(dstDir, 'VE11'), gt_string);

%CopyGadgetronInstallers('D:\installer\VE11C_gt', fullfile(dstDir, 'VE11C_gt'), chroot_string);
%CopyGadgetronInstallers('D:\installer\VE11C', fullfile(dstDir, 'VE11C'), chroot_string);
%CopyGadgetronInstallers('D:\installer\VE11S', fullfile(dstDir, 'VE11S'), chroot_string);

% CopyGadgetronInstallers('D:\installer\VE11B', fullfile(dstDir, 'VE11B'), chroot_string);
%CopyGadgetronInstallers('D:\installer\VE11P', fullfile(dstDir, 'VE11P'), chroot_string);

% CopyGadgetronInstallers('D:\installer\NXVA11A', fullfile(dstDir, 'NXVA11A'), chroot_string);
CopyGadgetronInstallers('D:\installer\NXVA20A', fullfile(dstDir, 'NXVA20A'), chroot_string);
CopyGadgetronInstallers('D:\installer\VE11E', fullfile(dstDir, 'VE11E'), chroot_string);
CopyGadgetronInstallers('D:\installer\gadgetron_fire\VE11E', fullfile(dstDir, 'VE11E'), chroot_string);


installFile = 'D:\Temp\install_chroot_image.sh';
if(isFileExist(installFile))
    copyfile('D:\Temp\install_chroot_image.sh', dstDir);
end

tarFile = fullfile('D:\Temp', [chroot_string '.tar.gz']);
if(isFileExist(tarFile))
    copyfile(tarFile, dstDir);
end

imgFile = fullfile('E:\temp', [chroot_string '.img']);
if(isFileExist(imgFile))
    copyfile(imgFile, dstDir);
end

winopen(dstDir);
