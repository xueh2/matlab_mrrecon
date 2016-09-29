function [mean_v, max_v, min_v, std_v] = ComputeMoCoPerformance_Perfusion_KLT_evaluation_mean(v)
% Compute mean/max/min/std from v

N = size(v, 1);
M = size(v, 2);

ind = find( abs(v(:))>0);

% s = [];
% for n=1:N
%     for m=n+1:M
%         if( abs(v(n,m))>0)
%             s = [s; v(n,m)];
%         end
%     end
% end

% s = abs(s);

if(isempty(ind))
    mean_v = 0;
    max_v = 0;
    min_v = 0;
    std_v = 0;
    return;
end

s = abs(v(ind(:)));

mean_v = mean(s(:));
max_v = max(s(:));
min_v = min(s(:));
std_v = std(s(:));
