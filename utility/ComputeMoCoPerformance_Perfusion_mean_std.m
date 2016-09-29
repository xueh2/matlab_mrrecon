function [m_v, std_v, m_v_KLT, std_v_KLT, m_v_progressive, std_v_progressive, m_v_direct, std_v_direct, pvalue_v_vKLT, pvalue_v_vPro, pvalue_vKLT_vPro, pvalue_vKLT_vDir] = ComputeMoCoPerformance_Perfusion_mean_std(v, v_KLT, v_progressive, v_direct, row)
% get the first pass results

m_v = squeeze(mean(v, 2));
std_v = squeeze(std(v, [], 2));

m_v_KLT = squeeze(mean(v_KLT, 2));
std_v_KLT = squeeze(std(v_KLT, [], 2));

m_v_progressive = squeeze(mean(v_progressive, 2));
std_v_progressive = squeeze(std(v_progressive, [], 2));

m_v_direct = squeeze(mean(v_direct, 2));
std_v_direct = squeeze(std(v_direct, [], 2));

[h, pvalue_v_vKLT] = ttest(v(row,:,:), v_KLT(row,:,:));
[h, pvalue_v_vPro] = ttest(v(row,:,:), v_progressive(row,:,:));
[h, pvalue_vKLT_vPro] = ttest(v_KLT(row,:,:), v_progressive(row,:,:));
[h, pvalue_vKLT_vDir] = ttest(v_KLT(row,:,:), v_direct(row,:,:));