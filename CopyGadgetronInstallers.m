function CopyGadgetronInstallers(compilingDir, dstDir, chroot_string)
% copy gadgetron installer to the dst folder
% CopyGadgetronInstallers(compilingDir, dstDir, chroot_string)

if (~isdir(dstDir))
    mkdir(dstDir);
end

[names, num] = findFILE(compilingDir, ['*' chroot_string '*.exe']);

if(num>0)
    for ii=1:num
        names{ii}
        [pathstr, name, ext] = fileparts(names{ii});
        copyfile(names{ii}, dstDir, 'f');
    end
end
