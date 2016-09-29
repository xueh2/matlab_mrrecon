function D = samplingPatternSparseMatrixRepresentation(kspace, idx_acq)
% D is a sparse sampling matrix sampledPoints = D*kspace(:) 
% idx_acq : the 1D indexes to mark the sampled points

Nfe = size(kspace, 1);
Npe = size(kspace, 2);
nCoil = size(kspace, 3);
N = numel(idx_acq);

dSize = [N Nfe*Npe*nCoil];
D = spalloc(dSize(1), dSize(2), N);
for n=1:N
    D(n, idx_acq(n)) = 1;
end