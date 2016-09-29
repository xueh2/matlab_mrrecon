function [mean_v, max_v, min_v, std_v] = ComputeMoCoPerformance_Perfusion_KLT_evaluation_mean_keyFrame(v, keyFrame, startFP, endFP)
% Compute mean/max/min/std from v

N = size(v, 1);

mean_v = zeros(N-1, 1);
max_v = zeros(N-1, 1);
min_v = zeros(N-1, 1);
std_v = zeros(N-1, 1);

for n=1:N-1
    p = v{n+1, 1};
    if(size(p, 3)>1)
        [mean_v(n, 1), max_v(n), min_v(n), std_v(n)] = ComputeMoCoPerformance_Perfusion_KLT_evaluation_mean_One(keyFrame, p(:,:,1), v{1}, startFP, endFP);
    else
        [mean_v(n, 1), max_v(n), min_v(n), std_v(n)] = ComputeMoCoPerformance_Perfusion_KLT_evaluation_mean_One(keyFrame, p, v{1}, startFP, endFP);
    end
end