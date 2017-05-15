
function timeUsed = CopyGadgetronDebugOutputOnRemote(host, debug_folder, dst_folder, deleteRemote)
%% copy gadgetron debug output
% timeUsed = CopyGadgetronDebugOutputOnRemote(host, debug_folder, dst_folder)

if(nargin < 4)
    deleteRemote = 0;
end

tic;

[key, user] = sshKeyLookup(host);

if( ~isFileExist(dst_folder) )
    mkdir(dst_folder);
end

% command = ['ssh -i ' key ' ' user '@' host ' "' 'rm -rf /home/' user '/Debug/DebugOutput.tar.gz' '"'];
% command
% dos(command, '-echo');
% 
% gt_command = ['tar -zcf /home/' user '/Debug/DebugOutput.tar.gz /home/' user '/Debug/DebugOutput/'];
% command = ['ssh -i ' key ' ' user '@' host ' "' gt_command '"'];
% command
% dos(command, '-echo');

command = ['pscp -r -i ' key '.ppk ' user '@' host ':' debug_folder ' ' dst_folder];
% command = ['scp -o StrictHostKeyChecking=no  -i ' key ' ' user '@' host ':' debug_folder ' ' dst_folder];
command
dos(command, '-echo');

% command = ['pscp -r -i ' key '.ppk ' user '@' host ':' '/home/' user '/Debug/DebugOutput.tar.gz ' dst_folder];
% command
% dos(command, '-echo');
% 
% gt_command = ['rm -rf ' '/home/' user '/Debug/DebugOutput.tar.gz '];
% command = ['ssh -i ' key ' ' user '@' host ' "' gt_command '"'];
% command
% dos(command, '-echo');
    
if(deleteRemote)
    gt_command = ['rm -rf ' debug_folder '/*.*'];
    command = ['ssh -i ' key ' ' user '@' host ' "' gt_command '"'];
    command
    dos(command, '-echo');    
end

% command = ['7z.exe e ' fullfile(dst_folder, 'DebugOutput.tar.gz') ];
% command
% dos(command, '-echo');
% 
% command = ['7z.exe e ' fullfile(dst_folder, 'DebugOutput.tar') ' -oDebugOutput' ];
% command
% dos(command, '-echo');
% 
% dos(['del /Q /F ' fullfile(dst_folder, 'DebugOutput.tar.gz')], '-echo');
% dos(['del /Q /F ' fullfile(dst_folder, 'DebugOutput.tar')], '-echo');

timeUsed = toc;