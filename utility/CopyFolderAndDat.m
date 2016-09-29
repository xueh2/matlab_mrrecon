function CopyFolderAndDat(homeBase, dstRoot)
% CopyFolderAndDat(home, dstRoot)
% copy all dat files with folder structure

[subdir1st, num1st] = FindSubDirs(homeBase);
dstDir = fullfile(dstRoot)
mkdir(dstDir);

[names, numF] = findFILE(homeBase, '*.dat');
if(numF>0)
    
    for kk=1:numF
        srcData = names{kk};
        
        [pathstr, dname, ext] = fileparts(srcData);
        dstData = fullfile(dstRoot, [dname ext]);
        
        doCopy = 1;
        if(isFileExist(dstData))
            a = dir(srcData);
            b = dir(dstData);
            if(a.bytes==b.bytes)
                doCopy = 0;
            end
        else
            doCopy = 1;
        end
        
        if(doCopy)
            disp(['Copy from ' srcData ' to ' dstData ' ... ']);
            copyfile(srcData, dstDir);        
        else
            disp(['Already exists - ' dstData ' ... ']);
        end
    end
end

srcDir = homeBase;
dstDir = dstRoot;

for ii=1:num1st
    if( ~isempty(strfind(subdir1st{ii}, 'res')) | ~isempty(strfind(subdir1st{ii}, 'grappa')) | ~isempty(strfind(subdir1st{ii}, 'generic')) | ~isempty(strfind(subdir1st{ii}, 'ref')) | ~isempty(strfind(subdir1st{ii}, 'reduced')) )
        continue;
    end
    currSrcDir = fullfile(srcDir, subdir1st{ii});
    currDstDir = fullfile(dstDir, subdir1st{ii});
    mkdir(currDstDir);              
    CopyFolderAndDat(currSrcDir, currDstDir);
end        

% for p=1:num1st
%     
%     home = fullfile(homeBase, subdir1st{p});
%     mkdir(fullfile(dstRoot, subdir1st{p}));
% 
%     [subdir2st, num2st] = FindSubDirs(home);
%     
%     [names, numF] = findFILE(home, '*.dat');
%     if(numF>0)
%         
%         dstDir = fullfile(dstRoot, src, folderName, data)
%         mkdir(dstDir);
%         
%         copyfile(fullfile(home, ['*.dat']), dstDir);        
%     end
%         
%     exDir = {'ICERecon', 'grappa', 'TXMapping', 'moco', 'mocoSyn', 'mocoPS' , 'mocoPSSyn', 'seg', 'magFitting', 'slep_res', 'slep_cloud_res', 'grappa_res', 'ICE', 'DebugOutput', 'PhantomT1maps', 'slep_flow_res', 'grappa_flow_res', 'slep_cloud_flow_res', 'slep_cloud_flow_res2'};
%     [subdir, num] = FindAllEndDirectoryExclusive(home, exDir)
% 
%     inds = strfind(home, '\');
%     src = home(inds(end)+1:end)
% 
%     for ii=1:num  
%         inds = strfind(subdir{ii}, '\');
%         v = subdir{ii};
%         if(~isempty(inds))
%             folderName = v(1:inds(end)-1)
%         else
%             folderName = []
%         end
% 
%         if(~isempty(inds))
%             data = v(inds(end)+1:end)   
%         else
%             data = v   
%         end
% 
%         dstDir = fullfile(dstRoot, src, folderName, data)
%         mkdir(dstDir);
% 
%         [names, numF] = findFILE(fullfile(home, folderName, data), '*.dat');
% 
%         if(numF>0)
%             copyfile(fullfile(home, folderName, data, ['*.dat']), dstDir);        
%         end
%     end
% end
