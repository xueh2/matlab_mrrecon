
function timeUsed = BuildGadgetronDockerPackageOnRemote(host, res_dir, docker_image, target_dir, create_installers)
%% build docker image on remote compute and copy them
% timeUsed = BuildGadgetronDockerPackageOnRemote('barbados', '/home/GADGETRON_CUDA75', 'xueh2/gadgetron_ubuntu1404_cuda_gtprep', '\\macau.nhlbi.nih.gov\data1\gadgetron_installation\20160405', 1)
% timeUsed = BuildGadgetronDockerPackageOnRemote('barbados', '/home/GADGETRON_CUDA55', 'hansenms/gadgetron_ubuntu1404_cuda_gtprep55', '\\137.187.134.135\share\Installer\20160512', 1)
% timeUsed = BuildGadgetronDockerPackageOnRemote('barbados', '/home/GADGETRON_CUDA75', 'xueh2/gadgetron_ubuntu1404_cuda_gtprep , '\\137.187.134.135\share\Installer\20160720', 1)
% timeUsed = BuildGadgetronDockerPackageOnRemote('denmark', '/home/GADGETRON_CUDA75', 'xueh2/gadgetron_ubuntu1604_cuda_gtprep', '\\137.187.134.135\share\Installer\20170525', 1)

if nargin < 1
    host = 'barbados'
    res_dir = '/home/GADGETRON_CUDA55'
    % docker_image = 'xueh2/gadgetron_ubuntu1404_cuda_gtprep'
    docker_image = 'hansenms/gadgetron_ubuntu1404_cuda_gtprep55'
    target_dir = '/e/temp'
    create_installers = 0
end

tic;

mkdir(target_dir);

[key, user] = sshKeyLookup(host);

gt_command = ['docker pull ' docker_image];

disp('pull latest image')
command = ['ssh -i ' key ' ' user '@' host ' "' gt_command '"'];
command
dos(command, '-echo');

% curr_t = datenum(clock);
% new_docker_image = [docker_image '_' num2str(curr_t)]
% disp('strip image')
% gt_command = ['cd ' res_dir ' && ./strip_docker_image ' docker_image ' ' new_docker_image];
% command = ['ssh -i ' key ' ' user '@' host ' "' gt_command '"'];
% command
% dos(command, '-echo');

disp('build package')
%gt_command = ['cd ' res_dir ' && sudo ./create_chroot_from_image ' new_docker_image ' 2048'];
gt_command = ['cd ' res_dir ' && sudo ./create_chroot_from_image ' docker_image ' 6144'];
command = ['ssh -i ' key ' ' user '@' host ' "' gt_command '"'];
command
dos(command, '-echo');

% disp('remove image')
% gt_command = ['cd ' res_dir ' && docker rmi ' new_docker_image];
% command = ['ssh -i ' key ' ' user '@' host ' "' gt_command '"'];
% command
% dos(command, '-echo');

disp('copy package')
gt_command = ['cd ' res_dir ' && ls -t | head -n1 '];
command = ['ssh -i ' key ' ' user '@' host ' "' gt_command '"'];
command
[s, fw] = dos(command, '-echo');
fw

[pathstr, filename, ext] = fileparts(fw);
gt_command = ['cd ' res_dir ' && ls -t | head -n1 '];
command = ['scp -i ' getenv('GADGETRON_KEY_FOLDER') '/gtuser_denmark_private ' user '@' host ':' res_dir '/' filename '.img ' target_dir];
command
dos(command, '-echo');

if(create_installers)
    command = [getenv('GADGETRON_SCRIPTS_FOLDER') '\compile_gadgetron_package_Sites_All_Versions.bat ' filename ' 3 14 1'];
    command
    dos(command, '-echo');
    
    disp(['target_dir is : ' target_dir]);
    disp(['filename is : ' filename]);
    CopyGadgetronInstallersAll(target_dir, filename, '3.14.1');
end

timeUsed = toc;
