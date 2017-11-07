
function CreateDirs(home, nameDir)
% create a nameDir directory in all subdirs in home

[subdirs, num] = FindAllDirectory(home);

for i = 1:num
    currentDir = [home '\' subdirs{i}];
    [success,meg,megID] = mkdir(currentDir, nameDir);
end

return;