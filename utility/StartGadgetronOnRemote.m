
function timeUsed = StartGadgetronOnRemote(host, port)
%% start gadgetron on remote server
% timeUsed = StartGadgetronOnRemote(host, port)

if nargin < 2
    port=9002;
end

tic;

[key, user] = sshKeyLookup(host);

gt_command = ['rm -rf /home/' user '/Debug/record_' num2str(port) '.txt'];
command = ['ssh -i ' key ' ' user '@' host ' "' gt_command '"'];
command
dos(command, '-echo');
    
gt_command = ['bash -c ''nohup /home/' user '/gt_scanner_setup_scripts/StartGadgetron ' num2str(port) ' > /home/' user '/Debug/record_' num2str(port) '.txt 2>&1 < /dev/null &'' '];
command = ['ssh -i ' key ' ' user '@' host ' "' gt_command '"'];
command
dos(command, '-echo');

timeUsed = toc;