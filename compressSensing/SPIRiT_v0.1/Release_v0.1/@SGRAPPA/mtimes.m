function res = mtimes(a,x)
% This method performs the SGRAPPA operator

idxacq = a.idxacq;
kernelIm = a.kernelIm;
nCoils = size(kernelIm,4);

res = zeros(size(x));
if a.adjoint
    xx = ifft2c(x);
    for n=1:nCoils
        tmpk = squeeze(conj(kernelIm(:,:,n,:)));
        res(:,:,n) = sum(tmpk.*xx,3); 
    end
    res = fft2c(res);
    res(idxacq) = x(idxacq);
else
    xx = ifft2c(x);
    for n=1:nCoils
        tmpk = kernelIm(:,:,:,n);
        res(:,:,n) = sum(tmpk.*xx,3);
    end
    res = fft2c(res);
    res(idxacq) = x(idxacq);
end
    
