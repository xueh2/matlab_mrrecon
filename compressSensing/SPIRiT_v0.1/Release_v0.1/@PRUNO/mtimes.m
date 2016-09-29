function res = mtimes(a,x)
% This method performs the PRUNO N'N*x operator
% (N'N)' = N'N, so for PRUNO kernel, the adjoint and original are same

kernelComposite = a.kernelComposite;
nCoils = size(kernelComposite,4);

res = zeros(size(x));
% if a.adjoint
% 
%     xx = ifft2c(x);
%     for n=1:nCoils
%         tmpk = squeeze(conj(kernelComposite(:,:,n,:)));
%         res(:,:,n) = sum(tmpk.*xx,3); 
%     end
%     res = fft2c(res);
% else
    xx = ifft2c(x);
    for n=1:nCoils
        tmpk = kernelComposite(:,:,:,n);
        res(:,:,n) = sum(tmpk.*xx,3);
    end
    res = fft2c(res);
% end
    
