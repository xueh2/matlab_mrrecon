function [mean_v, max_v, min_v, std_v] = ComputeMoCoPerformance_Perfusion_KLT_evaluation_mean_all(v)
% Compute mean/max/min/std from v

N = size(v, 1);
M = size(v, 2);

mean_v = zeros(N, M);
max_v = zeros(N, M);
min_v = zeros(N, M);
std_v = zeros(N, M);

for n=1:N
    for m=1:M
        p = v{n, m};
        if(size(p, 3)>1)
            [mean_v(n, m), max_v(n, m), min_v(n, m), std_v(n, m)] = ComputeMoCoPerformance_Perfusion_KLT_evaluation_mean(p(:,:,1));
        else
            [mean_v(n, m), max_v(n, m), min_v(n, m), std_v(n, m)] = ComputeMoCoPerformance_Perfusion_KLT_evaluation_mean(p);
        end
    end
end