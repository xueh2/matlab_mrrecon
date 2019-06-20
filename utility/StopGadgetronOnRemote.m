
function timeUsed = StopGadgetronOnRemote(host, port)
%% stop gadgetron on remote server
% timeUsed = StopGadgetronOnRemote(host, port)

tic;
[key, user] = sshKeyLookup(host);

is_remote_computer = IsRemoteComputer(host);

if(~is_remote_computer)
    gt_command = ['/home/' user '/gt_scanner_setup_scripts/StopGadgetronService >> /home/' user '/Debug/record_' num2str(port) '.txt 2>&1 < /dev/null &'];
    dos(gt_command, '-echo');
else

    gt_command = ['sh -c ''/home/' user '/gt_scanner_setup_scripts/StopGadgetronService >> /home/' user '/Debug/record_' num2str(port) '.txt 2>&1 < /dev/null &'' '];

    % command = ['ssh -i ' key ' ' user '@' host ' "' gt_command '"']

    is_remote_computer = IsRemoteComputer(host);   
    if(~is_remote_computer)
        command = gt_command;
    else
        command = ['ssh ' user '@' host ' "' gt_command '"']
    end
    dos(command, '-echo');
end

timeUsed = toc;