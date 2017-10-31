
function [target, source, rah, numofh] = ParseDofName(dofname)

target = [];
source = [];
rah = [];
numofh = 0;

index = find(dofname == '-');
len = length(index);

target = dofname(1:index(1)-1);
source = dofname( (index(1)+1):(index(2)-1));
rah = dofname( (index(2)+1):(index(2)+4));
if ( len > 2 )
    numofh = num2str( dofname(index(3)+1) );
end

index = strfind(target,'_t2');
target = target(1:index-1);

index = strfind(source,'_t2');
source = source(1:index-1);

return;