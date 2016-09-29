function kernel = spiritCalibration(kCalib, kSize, thresReg)

[sx,sy,nCoil] = size(kCalib);

A = [];
D = [];

dummyK = zeros(kSize(1),kSize(2)); 
dummyK((end+1)/2,(end+1)/2) = 1;
idxCentre = find(dummyK==1);
idxUsed = find(dummyK==0);

ind = [1:idxCentre-1 idxCentre+1:kSize(1)*kSize(2)];

% A*kernel = D

% composite the data matrix
for n=1:nCoil
	tmp  = im2col(kCalib(:,:,n),kSize,'sliding').';    
	A = [A, tmp(:, ind(:))];
    D = [D tmp(:, idxCentre)];
end

% if ( hasGPU() )
%     tic
%     A = gpuArray(single(A));
%     D = gpuArray(single(D));
%     toc
% 
%     tic
%     AtA = A'*A;
%     Aty = A'*D;
%     toc
% 
%     tic
%     dTraceR = trace(real(AtA));
%     nRow = size(AtA, 1);
%     lamda  = thresReg * dTraceR / nRow;
%     for n=1:nRow
%         AtA(n,n) = real(AtA(n,n)) + lamda;
%     end
%     toc
%     
%     tic
%     rawkernel = AtA\Aty;
%     toc
%     
%     rawkernel = gather(rawkernel);
% else
    AtA = A'*A;
    Aty = A'*D;
    
    dTraceR = trace(real(AtA));
    nRow = size(AtA, 1);
    lamda  = thresReg * dTraceR / nRow;
    for n=1:nRow
        AtA(n,n) = real(AtA(n,n)) + lamda;
    end
    
    rawkernel = AtA\Aty;

% end

dummyK = zeros(kSize(1),kSize(2),nCoil); 
dummyK((end+1)/2,(end+1)/2,:) = 1;
idxUsed = find(dummyK==0);

aKernel = zeros(kSize(1),kSize(2),nCoil);
kernel = zeros(kSize(1),kSize(2),nCoil,nCoil);
for n=1:nCoil
    aKernel(idxUsed) = rawkernel(:,n);
    kernel(:,:,:,n) = aKernel;
end
