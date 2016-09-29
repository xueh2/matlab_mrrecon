function kernel = spiritCalibration3D(kCalib, kSize, thresReg)
% kCalib: [COL LIN PAR CHA]
% ksize: [kfe kpe kpar]
% kernel: [kfe kpe kpar srcCHA dstCHA]
% kernel is arranged as a 5D array
% This calibration is slightly different from Lustig's spirit
% the central point of every kSize patch across all coils does not
% contribute to data fitting

[sx,sy,sz,nCoil] = size(kCalib);

A = [];
D = [];

dummyK = zeros(kSize(1),kSize(2),kSize(3),nCoil); 
dummyK((end+1)/2,(end+1)/2,(end+1)/2,:) = 1;
idxCentre = find(dummyK(:)==1);
idxUsed = find(dummyK(:)==0);

% A*kernel = D

% composite the data matrix
startFE = ceil(kSize(1)/2);
endFE = sx - floor(kSize(1)/2);

startPE = ceil(kSize(2)/2);
endPE = sy - floor(kSize(2)/2);

startPar = ceil(kSize(3)/2);
endPar = sz - floor(kSize(3)/2);

M = (endFE-startFE+1)*(endPE-startPE+1)*(endPar-startPar+1);
N = (kSize(1)*kSize(2)*kSize(3)-1)*nCoil;
K = nCoil;

A = zeros(M, N);
D = zeros(M, K);

halfKSize = floor(kSize/2);

index = 1;
for par=startPar:endPar
    for pe=startPE:endPE
        for fe=startFE:endFE            
            srcData = kCalib(fe-halfKSize(1):fe+halfKSize(1), pe-halfKSize(2):pe+halfKSize(2), par-halfKSize(3):par+halfKSize(3),:);            
            A(index,:) = srcData(idxUsed(:));            
            D(index, :) = srcData(idxCentre(:));
            index = index + 1;
        end
    end
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

dummyK = zeros(kSize(1),kSize(2),kSize(3),nCoil); 
dummyK((end+1)/2,(end+1)/2,(end+1)/2,:) = 1;
idxUsed = find(dummyK==0);

aKernel = zeros(kSize(1),kSize(2),kSize(3),nCoil);
kernel = zeros(kSize(1),kSize(2),kSize(3),nCoil,nCoil);
for n=1:nCoil
    aKernel(idxUsed) = rawkernel(:,n);
    kernel(:,:,:,:,n) = aKernel;
end
