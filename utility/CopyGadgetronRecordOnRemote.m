
function timeUsed = CopyGadgetronRecordOnRemote(host, port, dstFile)
%% copy gadgetron record to local
% timeUsed = CopyGadgetronRecordOnRemote(host, port, dstFile)

tic;

[key, user] = sshKeyLookup(host);

if(strcmp(host, 'localhost')==1)
    command = ['cp ' '/home/' user '/Debug/record_' num2str(port) '.txt ' dstFile  ];
else  
    command = ['scp ' user '@' host ':' '/home/' user '/Debug/record_' num2str(port) '.txt ' dstFile  ];
end

dos(command, '-echo');

timeUsed = toc;