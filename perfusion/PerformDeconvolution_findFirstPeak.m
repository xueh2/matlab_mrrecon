function [F, ind] = PerformDeconvolution_findFirstPeak(impulse, maxDistRatio)
% find first peak of impulse response

%% set up common parameters
for n=1:numel(impulse)-1
    v = impulse(n);
    if ( v>0 && impulse(n)>impulse(n+1) )
        break;
    end                               
end

F = -1;
ind = n;
if ( n < numel(impulse)*maxDistRatio && v>0 )
    F = v;
end

