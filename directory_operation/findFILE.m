
function [names, num] = findFILE(home, suffix)
% find all files with the suffix in directory home
% [names, num] = findFILE(home, suffix)

if (exist(home) ~= 7 )
    names = [];
    num = 0;
    return;
end

currDir=pwd;
cd(home)
pp = dir(suffix);
num = numel(pp);

names = cell(num, 1);

for i=1:num
    names{i} = fullfile(home, pp(i).name);
end
cd(currDir)