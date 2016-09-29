function [kernel, rawkernel] = PFConvKernelCalibration(kSrc, kDst, kSize, thresReg)

[sx,sy,srcCHA] = size(kSrc);
[sx,sy,dstCHA] = size(kDst);

A = [];
D = [];

dummyK = zeros(kSize(1),kSize(2)); 
idxUsed = find(dummyK==0);

dummyK((end+1)/2,(end+1)/2) = 1;
idxCentre = find(dummyK==1); 

ind = [1:kSize(1)*kSize(2)];

% A*kernel = D

% composite the data matrix
for n=1:srcCHA
	tmp  = im2col(kSrc(:,:,n),kSize,'sliding').';    
	A = [A, tmp(:, ind(:))];
end

for n=1:dstCHA
    tmp2  = im2col(kDst(:,:,n),kSize,'sliding').';    
    D = [D tmp2(:, idxCentre)];
end

AtA = A'*A;
Aty = A'*D;

dTraceR = trace(real(AtA));
nRow = size(AtA, 1);
lamda  = thresReg * dTraceR / nRow;
for n=1:nRow
    AtA(n,n) = real(AtA(n,n)) + lamda;
end

rawkernel = AtA\Aty;

aKernel = zeros(kSize(1),kSize(2),srcCHA);
kernel = zeros(kSize(1),kSize(2),srcCHA,dstCHA);
for n=1:dstCHA
    aKernel = rawkernel(:,n);
    kernel(:,:,:,n) = reshape(aKernel, [kSize(1),kSize(2),srcCHA]);
end
