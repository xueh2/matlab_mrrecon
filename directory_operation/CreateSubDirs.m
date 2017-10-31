
function [subdirs, num] = CreateSubDirs(home, subdirs)
% create a nameDir directory in all subdirs in home

num = numel(subdirs);

for i=1:num
    mkdir(fullfile(home, subdirs{i}));
end