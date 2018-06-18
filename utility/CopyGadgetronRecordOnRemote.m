
function timeUsed = CopyGadgetronRecordOnRemote(host, port, dstFile)
%% copy gadgetron record to local
% timeUsed = CopyGadgetronRecordOnRemote(host, port, dstFile)

tic;

[key, user] = sshKeyLookup(host);

% command = ['pscp -i ' key '.ppk ' user '@' host ':' '/home/' user '/Debug/record_' num2str(port) '.txt ' dstFile  ];
command = ['scp ' user '@' host ':' '/home/' user '/Debug/record_' num2str(port) '.txt ' dstFile  ];
dos(command, '-echo');

timeUsed = toc;