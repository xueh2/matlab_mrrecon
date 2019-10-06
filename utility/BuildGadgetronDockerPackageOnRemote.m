
function timeUsed = BuildGadgetronDockerPackageOnRemote(host, res_dir, docker_image, target_dir, create_installers, remove_gtprep_xml, strip_image, chroot_size)
%% build docker image on remote compute and copy them
% create_installers = 1
% remove_gtprep_xml = 0
% strip_image = 0
% chroot_size = round(3.75*1024)
% timeUsed = BuildGadgetronDockerPackageOnRemote('grenada', '/home/GADGETRON', 'gadgetronnhlbi/gadgetron_ubuntu1604_gtprep_stripped', '\\hl-share\RawMRI\Lab-Kellman\Share\Installers\20190420', create_installers, remove_gtprep_xml, strip_image, chroot_size)
%
% create_installers = 1
% remove_gtprep_xml = 0
% strip_image = 0
% chroot_size = 3*1024
% timeUsed = BuildGadgetronDockerPackageOnRemote('grenada', '/home/GADGETRON', 'gadgetron/ubuntu_1804_no_cuda', '\\hl-share\RawMRI\Lab-Kellman\Share\Installers\20190420', create_installers, remove_gtprep_xml, strip_image, chroot_size)

if nargin < 1
    host = 'denmark'
    res_dir = '/home/GADGETRON_CUDA80'
    docker_image = 'xueh2/gadgetron_ubuntu1404_cuda_gtprep'
%     docker_image = 'xueh3/gadgetron_ubuntu1604_gtprep'
    target_dir = '/e/temp'
    create_installers = 0
end

if nargin < 6
    remove_gtprep_xml = 0
end

if nargin < 7
    strip_image = 0;
end

if nargin < 8
    chroot_size = 4.0*1024;
end

tic;

mkdir(target_dir);

[key, user] = sshKeyLookup(host);
user = 'xueh2'

gt_command = ['docker pull ' docker_image];

disp('pull latest image')
command = ['ssh ' user '@' host ' "' gt_command '"'];
command
dos(command, '-echo');

if( strip_image )
    curr_t = datenum(clock);
    new_docker_image = [docker_image '_' num2str(curr_t)]
    disp('strip image')
    gt_command = ['cd ' res_dir ' && ./strip_docker_image ' docker_image ' ' new_docker_image];
    command = ['ssh '  user '@' host ' "' gt_command '"'];
    command
    dos(command, '-echo');
else
    new_docker_image = docker_image;
end

disp('build package')
%gt_command = ['cd ' res_dir ' && sudo ./create_chroot_from_image ' new_docker_image ' 2048'];
if(remove_gtprep_xml)
    % gt_command = ['cd ' res_dir ' && sudo ./create_chroot_from_image ' docker_image ' 6144 GTPrep_*D*.xml'];
    gt_command = ['cd ' res_dir ' && sudo ./create_chroot_from_image ' new_docker_image ' ' num2str(chroot_size) ' GTPrep_*DT*.xml'];
else
    % gt_command = ['cd ' res_dir ' && sudo ./create_chroot_from_image ' docker_image ' 6144'];
    gt_command = ['cd ' res_dir ' && sudo ./create_chroot_from_image ' new_docker_image ' ' num2str(chroot_size) ' '];
end

command = ['ssh '  user '@' host ' "' gt_command '"'];
command
dos(command, '-echo');

if( strip_image )
    disp('remove stripped image')
    gt_command = ['cd ' res_dir ' && docker rmi ' new_docker_image];
    command = ['ssh '  user '@' host ' "' gt_command '"'];
    command
    dos(command, '-echo');
end

disp('copy package')
gt_command = ['cd ' res_dir ' && ls -t | head -n1 '];
command = ['ssh '  user '@' host ' "' gt_command '"'];
command
[s, fw] = dos(command, '-echo');

ind = strfind(fw, 'gadgetron');
fw = fw(ind:end);
ind = strfind(fw, '.img');
fw = fw(1:ind+3);
[pathstr, filename, ext] = fileparts(fw);

gt_command = ['cd ' res_dir ' && zip ' filename '.img.zip ' filename '.img' ];
command = ['ssh '  user '@' host ' "' gt_command '"'];
command
tic; [s, fw] = dos(command, '-echo'); toc

command = ['scp ' user '@' host ':' res_dir '/' filename '.img.zip ' target_dir];
command
tic; dos(command, '-echo'); toc

command = ['scp ' user '@' host ':' res_dir '/' filename '.img ' target_dir];
command
tic; dos(command, '-echo'); toc

if(create_installers)
        
    if(remove_gtprep_xml)
        PerformGadgetronRecon_Encrypt_GtPrep_config('D:\gtuser\mrprogs\install\share\gadgetron\config', 'D:\gadgetron_scanner_setup\VE11\config');
        PerformGadgetronRecon_Encrypt_GtPrep_config('D:\gtuser\mrprogs\install\share\gadgetron\config', 'D:\gadgetron_scanner_setup\VE11C\config');    
        PerformGadgetronRecon_Encrypt_GtPrep_config('D:\gtuser\mrprogs\install\share\gadgetron\config', 'D:\gadgetron_scanner_setup\VE11E\config');    

        command = [getenv('GADGETRON_SCRIPTS_FOLDER') '\compile_gadgetron_package_MARS_All_Versions.bat ' filename ' 3 17 0'];
    else
        delete('D:\gtuser\mrprogs\gt_scanner_setup\VE11\config\*.xml')
        delete('D:\gtuser\mrprogs\gt_scanner_setup\VE11C\config\*.xml')
        delete('D:\gtuser\mrprogs\gt_scanner_setup\VE11E\config\*.xml')
        
        command = [getenv('GADGETRON_SCRIPTS_FOLDER') '\compile_gadgetron_package_Sites_All_Versions.bat ' filename ' 3 17 0'];
    end
       
    command
    tic; dos(command, '-echo'); toc
    
    disp(['target_dir = ' target_dir]);
    disp(['filename = ' filename]);
    CopyGadgetronInstallersAll(target_dir, filename, '3.17.0');
end

timeUsed = toc;
