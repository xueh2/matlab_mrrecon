
function stopflag = StopInflation(ThresRatio, target_stats, source_stats)
% stop inflationg if all measures are close in quantity

num = length(target_stats);

errors = abs(source_stats - target_stats);
flags = (errors <= ThresRatio*target_stats) | (source_stats <= target_stats);

if ( isempty(find(flags==0)) == 1 )
    stopflag = 1;
else
    stopflag = 0;
end

return