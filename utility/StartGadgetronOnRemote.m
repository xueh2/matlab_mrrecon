
function timeUsed = StartGadgetronOnRemote(host, port)
%% start gadgetron on remote server
% timeUsed = StartGadgetronOnRemote(host, port)

if nargin < 2
    port=9002;
end

tic;

[key, user] = sshKeyLookup(host);

gt_command = ['rm -rf /home/' user '/Debug/record_' num2str(port) '.txt'];
if(isempty(strfind(host, 'localhost')))
    command = ['ssh -o StrictHostKeyChecking=no ' user '@' host ' "' gt_command '"'];
else
    command = gt_command;
end
command
dos(command, '-echo');
    
gt_command = ['bash -c ''nohup /home/' user '/gt_scanner_setup_scripts/StartGadgetron ' num2str(port) ' 137.187.134.184 ' ' > /home/' user '/Debug/record_' num2str(port) '.txt 2>&1 < /dev/null &'' '];

if(isempty(strfind(host, 'localhost')))
    command = ['ssh -o StrictHostKeyChecking=no ' user '@' host ' "' gt_command '"'];
else
    command = gt_command;
end
command
dos(command, '-echo');

timeUsed = toc;