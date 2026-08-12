
function findAndMoveMeasDat(home, ext_str)
% find and move meas data files

if(nargin<1)
    home = '.';
end

if(nargin<2)
    ext_str = '.dat';
end

[names, num] = findFILE(home, ['*' ext_str])

use_h5 = 0;
if(num==0)
    [names, num] = findFILE(home, '*.h5')
    use_h5 = 1;
end

for ii=1:num
   
    [path, name, ext] = fileparts(names{ii});
    
    dstDir = fullfile(home, name);
    mkdir(dstDir);
    
    if(use_h5)
%         dos(['mv ' names{ii} ' ' fullfile(home, name, [name '.h5'])]);
         movefile(names{ii}, fullfile(home, name, [name '.h5']));
    else
%         dos(['mv ' names{ii} ' ' fullfile(home, name, [name '.dat'])]);
         movefile(names{ii}, fullfile(home, name, [name ext_str]));
    end
end