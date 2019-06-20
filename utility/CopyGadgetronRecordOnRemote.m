
function timeUsed = CopyGadgetronRecordOnRemote(host, port, dstFile)
%% copy gadgetron record to local
% timeUsed = CopyGadgetronRecordOnRemote(host, port, dstFile)

tic;

[key, user] = sshKeyLookup(host);

is_remote_computer = IsRemoteComputer(host);

if(~is_remote_computer)
    command = ['cp ' '/home/' user '/Debug/record_' num2str(port) '.txt ' dstFile  ];
else  
    command = ['scp ' user '@' host ':' '/home/' user '/Debug/record_' num2str(port) '.txt ' dstFile  ];
end

dos(command, '-echo');

timeUsed = toc;