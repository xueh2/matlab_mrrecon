
function timeUsed = CopyGadgetronDebugOutputOnRemote(host, debug_folder, dst_folder, deleteRemote)
%% copy gadgetron debug output
% timeUsed = CopyGadgetronDebugOutputOnRemote(host, debug_folder, dst_folder)

if(nargin < 4)
    deleteRemote = 0;
end

tic;

[key, user] = sshKeyLookup(host);

% if( ~isFileExist(dst_folder) )
%     mkdir(dst_folder);
% end

% mkdir(fullfile(dst_folder, 'DebugOutput'));

gt_command = ['mkdir -p ' dst_folder];
command = ['ssh -o TCPKeepAlive=yes -o ServerAliveInterval=15 -o ServerAliveCountMax=3 ' user '@' host ' "' gt_command '"'];
command
dos(command, '-echo');
    
% command = ['ssh -i ' key ' ' user '@' host ' "' 'rm -rf /home/' user '/Debug/DebugOutput.tar.gz' '"'];
% command
% dos(command, '-echo');
% 
% gt_command = ['tar -zcf /home/' user '/Debug/DebugOutput.tar.gz /home/' user '/Debug/DebugOutput/'];
% command = ['ssh -i ' key ' ' user '@' host ' "' gt_command '"'];
% command
% dos(command, '-echo');

% if(isunix())
%     mkdir(fullfile(dst_folder, 'DebugOutput'));
%     command = ['scp -q -o StrictHostKeyChecking=no ' user '@' host ':' debug_folder '/* ' dst_folder '/DebugOutput/'];
%     command
    
    gt_command = ['cp -r ' debug_folder '/* ' dst_folder];
    command = ['ssh -o TCPKeepAlive=yes -o ServerAliveInterval=15 -o ServerAliveCountMax=3 ' user '@' host ' "' gt_command '"'];
    command
    tic; dos(command, '-echo'); copy_duration = toc;

%     command = ['scp -q -o StrictHostKeyChecking=no ' user '@' host ':' debug_folder '/incomingKSpace_REAL* ' dst_folder '/DebugOutput/'];
%     tic; dos(command, '-echo'); copy_duration = toc;
%     command = ['scp -q -o StrictHostKeyChecking=no ' user '@' host ':' debug_folder '/incomingKSpace_IMAG* ' dst_folder '/DebugOutput/'];
%     tic; dos(command, '-echo'); copy_duration = toc;
%     
%     command = ['scp -q -o StrictHostKeyChecking=no ' user '@' host ':' debug_folder '/coilMap_* ' dst_folder '/DebugOutput/'];
%     tic; dos(command, '-echo'); copy_duration = toc;
%     command = ['scp -q -o StrictHostKeyChecking=no ' user '@' host ':' debug_folder '/complexIm_* ' dst_folder '/DebugOutput/'];
%     tic; dos(command, '-echo'); copy_duration = toc;
%     command = ['scp -q -o StrictHostKeyChecking=no ' user '@' host ':' debug_folder '/CmrSpatioTemporalFilter_* ' dst_folder '/DebugOutput/'];
%     tic; dos(command, '-echo'); copy_duration = toc;
%     command = ['scp -q -o StrictHostKeyChecking=no ' user '@' host ':' debug_folder '/input_for_fil_* ' dst_folder '/DebugOutput/'];
%     tic; dos(command, '-echo'); copy_duration = toc;
%     command = ['scp -q -o StrictHostKeyChecking=no ' user '@' host ':' debug_folder '/res_after_fil_* ' dst_folder '/DebugOutput/'];
%     tic; dos(command, '-echo'); copy_duration = toc;
%     command = ['scp -q -o StrictHostKeyChecking=no ' user '@' host ':' debug_folder '/res_afterunwrapping_* ' dst_folder '/DebugOutput/'];
%     tic; dos(command, '-echo'); copy_duration = toc;
%     command = ['scp -q -o StrictHostKeyChecking=no ' user '@' host ':' debug_folder '/kspace_before_POCS_* ' dst_folder '/DebugOutput/'];
%     tic; dos(command, '-echo'); copy_duration = toc;
%     command = ['scp -q -o StrictHostKeyChecking=no ' user '@' host ':' debug_folder '/kspace_after_POCS_* ' dst_folder '/DebugOutput/'];
%     tic; dos(command, '-echo'); copy_duration = toc;
%     command = ['scp -q -o StrictHostKeyChecking=no ' user '@' host ':' debug_folder '/unwarppedIm_* ' dst_folder '/DebugOutput/'];
%     tic; dos(command, '-echo'); copy_duration = toc;
  
% else    
%     command = ['pscp -r ' user '@' host ':' debug_folder ' ' dst_folder];
% end

% command
% tic; dos(command, '-echo'); copy_duration = toc;
disp(['Copy debug output ' num2str(copy_duration)]);

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
    command = ['ssh -o TCPKeepAlive=yes -o ServerAliveInterval=15 -o ServerAliveCountMax=3 ' user '@' host ' "' gt_command '"'];
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