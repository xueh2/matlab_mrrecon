
function Resample_TargetGridRun_Binary(target, source, output)
% resample the source in the target grid

command = ['resample_TargetGrid' ' ' target ' ' source ' ' output];
command
[s, w] = dos(command, '-echo');

return;