
function timeUsed = UpdateGadgetronOnRemote(host, clean, branch, gt_branch, gtprep_branch)
%% update gadgetron on remote server
% timeUsed = UpdateGadgetronOnRemote(host)
% UpdateGadgetronOnRemote('barbados')
% UpdateGadgetronOnRemote('denmark')
% UpdateGadgetronOnRemote('palau')
% UpdateGadgetronOnRemote('samoa')
% UpdateGadgetronOnRemote('andorra')

if(nargin<2)
    clean = 0;
end

if(nargin<3)
    branch = 'master';
end


if(nargin<4)
    gt_branch = 'master';
end

if(nargin<5)
    gtprep_branch = 'master';
end

tic;

[key, user] = sshKeyLookup(host);

  
gt_command = ['/home/' user '/gt_scanner_setup_scripts/UpdateGadgetron ' num2str(clean) ' ' branch ' ' gt_branch ' ' gtprep_branch];
command = ['ssh -i ' key ' ' user '@' host ' "' gt_command '"'];
command
dos(command, '-echo');

timeUsed = toc;