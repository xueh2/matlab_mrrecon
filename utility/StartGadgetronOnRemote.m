
function timeUsed = StartGadgetronOnRemote(host, port)
%% start gadgetron on remote server
% timeUsed = StartGadgetronOnRemote(host, port)

if nargin < 2
    port=9002;
end

tic;

[key, user] = sshKeyLookup(host);

is_remote_computer = IsRemoteComputer(host);

gt_command = ['rm -rf /home/' user '/Debug/record_' num2str(port) '.txt'];
if(is_remote_computer)
    command = ['ssh -o ConnectTimeout=10 -o StrictHostKeyChecking=no -o TCPKeepAlive=yes -o ServerAliveInterval=15 -o ServerAliveCountMax=3 ' user '@' host ' "' gt_command '"'];
else
    command = gt_command;
end
command
dos(command, '-echo');
  
if(is_remote_computer)
    gt_command = ['/home/' user '/gt_scanner_setup_scripts/StartGadgetron ' num2str(port) ' 8899 9988 137.187.135.157 ' ' > /home/' user '/Debug/record_' num2str(port) '.txt 2>&1 < /dev/null &'];    
else
    % gt_command = ['/home/' user '/gt_scanner_setup_scripts/StartGadgetron ' num2str(port) ' 8899 9988 137.187.135.157 > /home/' user '/Debug/record_' num2str(port) '.txt 2>&1 < /dev/null &'];
    % gt_command = ['nohup /home/' user '/gt_scanner_setup_scripts/StartGadgetron ' num2str(port) ' 8899 9988 137.187.135.157 ' ' > /home/' user '/Debug/record_' num2str(port) '.txt 2>&1 < /dev/null &'];    
    gt_command = ['gadgetron -p ' num2str(port) ' > /home/' user '/Debug/record_' num2str(port) '.txt 2>&1 < /dev/null &'];    
end

if(is_remote_computer)
    command = ['ssh -o ConnectTimeout=10 -o StrictHostKeyChecking=no -o TCPKeepAlive=yes -o ServerAliveInterval=15 -o ServerAliveCountMax=3 ' user '@' host ' "' gt_command '"'];
else
    command = gt_command;
end
command
dos(command, '-echo');

timeUsed = toc;