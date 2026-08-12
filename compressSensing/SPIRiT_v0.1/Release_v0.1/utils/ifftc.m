function res = ifftc(x,dim)
%res = ifftc(x,dim)
res = sqrt(size(x,dim))*ifftshift(ifft(fftshift(x,dim),[],dim),dim);

