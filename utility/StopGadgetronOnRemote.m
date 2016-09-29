
function timeUsed = StopGadgetronOnRemote(host, port)
%% stop gadgetron on remote server
% timeUsed = StopGadgetronOnRemote(host, port)

tic;

[key, user] = sshKeyLookup(host);

gt_command = ['sh -c ''/home/' user '/gt_scanner_setup_scripts/StopGadgetronService >> /home/' user '/Debug/record_' num2str(port) '.txt 2>&1 < /dev/null &'' '];

command = ['ssh -i ' key ' ' user '@' host ' "' gt_command '"']
dos(command, '-echo');

timeUsed = toc;