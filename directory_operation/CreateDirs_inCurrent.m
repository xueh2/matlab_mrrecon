
function CreateDirs_inCurrent(home, nameDir)
% create a nameDir directory in home

[success,meg,megID] = mkdir(home, nameDir);

return;