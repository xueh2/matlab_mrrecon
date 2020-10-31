
function Resample_TargetGridRun(target, source, output)
% resample the source in the target grid

command = ['resample_TargetGrid' ' ' target ' ' source ' ' output ' '  '-bspline' ];
command
[s, w] = dos(command, '-echo');

return;