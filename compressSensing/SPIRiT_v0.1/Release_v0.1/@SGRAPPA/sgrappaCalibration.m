function kernel = sgrappaCalibration(kCalib,kSize,acquiredLines,thresReg)
% compute the sgrappa kernel
% kspace: the fullkspace (Nfe Npe Coil), calibration area
% kSize : kernel size, kSize(1) is along row direction, kSize(2) along column
% Sampling pattern is assumed to be regularly sampled along PE
% acquiredLines : the index of originally acquired lines in the kSize(1), first line is line 1
% fitted line is the central line floor(kSize(2)/2) which should be the acquired lines as well
% NOTE: for the sgrappa point view, acquiredLines are actually lines to be estimated; it is acquired in the original kspace
% thresReg : threshold for regularization
% show : if 1, plot a figure to show the eigenvalue
% kernel : a kSize(1)*kSize(2)*nCoil*r matrix, corresponding to nCoil sgrappa kernels

disp(['calibration region size : ' num2str(size(kCalib))]);
disp(['kernel size : ' num2str(kSize)]);

[sx,sy,nCoil] = size(kCalib);

A = [];
D = [];

dummyK = zeros(kSize(1),kSize(2)); 
dummyK((end+1)/2,(end+1)/2) = 1;
idxCentre = find(dummyK==1);

for r=1:numel(acquiredLines)
    dummyK(:, acquiredLines(r)) = 1;
end
idxUsed = find(dummyK==0);
    
%  ind = [1:idxCentre-1 idxCentre+1:kSize(1)*kSize(2)];

% A*kernel = D

% composite the data matrix
for n=1:nCoil
	tmp  = im2col(kCalib(:,:,n),kSize,'sliding').';    
	A = [A, tmp(:, idxUsed(:))];
    D = [D tmp(:, idxCentre)];
end
AtA = A'*A;
Aty = A'*D;

% tic
% rawkernel = inv_tikhonov_IcePAT(AtA, thresReg)*Aty;
% toc

% tic
dTraceR = trace(real(AtA));
nRow = size(AtA, 1);
lamda  = thresReg * dTraceR / nRow;
for n=1:nRow
    AtA(n,n) = real(AtA(n,n)) + lamda;
end
rawkernel = AtA\Aty;
% toc

% norm(AtA*rawkernel-Aty)

dummyK = zeros(kSize(1),kSize(2),nCoil); 
for r=1:numel(acquiredLines)
    dummyK(:, acquiredLines(r), :) = 1;
end
idxUsed = find(dummyK==0);

aKernel = zeros(kSize(1),kSize(2),nCoil);
kernel = zeros(kSize(1),kSize(2),nCoil,nCoil);
for n=1:nCoil
    aKernel(idxUsed) = rawkernel(:,n);
    kernel(:,:,:,n) = aKernel;
end
