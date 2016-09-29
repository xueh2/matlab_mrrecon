function v_FP = ComputeMoCoPerformance_Perfusion_Get_FirstPass_Result(v, start_FP, end_FP)
% get the first pass results

N = numel(v);

v_FP = cell(N, 1);

num = end_FP-start_FP+1;
ind = find(v{1}>=start_FP & v{1}<=end_FP);

p = v{1};
v_FP{1} = p(ind(:), 1);

for kk=2:N   
    p = v{kk};
    v_FP{kk} = p(ind(:), ind(:), :);    
end