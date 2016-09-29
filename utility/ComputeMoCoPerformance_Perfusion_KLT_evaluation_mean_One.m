function [mean_v, max_v, min_v, std_v] = ComputeMoCoPerformance_Perfusion_KLT_evaluation_mean_One(keyFrame, v, frames, startFP, endFP)
% Compute mean/max/min/std from v

ind = find(frames(:)==keyFrame);
if(isempty(ind))
    error('frames do not have keyFrame');
end

if(startFP>0 && endFP>0)
    indS = find(frames(:)>=startFP);
    indS = indS(1);

    indE = find(frames(:)<=endFP);
    indE = indE(end);
end

indAll = [1:ind(1)-1 ind(1)+1:size(v, 2)];

if(startFP>0 && endFP>0)
    indAll = intersect(indAll, indS:indE);
end

vKey = v(ind(1), indAll(:));

s = abs(vKey(:));

mean_v = mean(s(:));
max_v = max(s(:));
min_v = min(s(:));
std_v = mean(s(:));