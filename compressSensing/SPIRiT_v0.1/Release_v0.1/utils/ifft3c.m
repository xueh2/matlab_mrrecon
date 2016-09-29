function res = ifft3c(x)

s = size(x);
fctr = size(x,1)*size(x,2)*size(x,3);

res = x;

if ( numel(s) == 3 ) % Nfe Npe Npar
    res = ifftshift(x);    
    res = ifftn(res);    
    res = fftshift(res);
    
    res = sqrt(fctr)*res;
end

if ( numel(s) == 4 ) % Nfe Npe Npar Cha
    for n=1:size(x,4)        
        res(:,:,:,n) = ifftshift(x(:,:,:,n));
        res(:,:,:,n) = ifftn(res(:,:,:,n));
        res(:,:,:,n) = sqrt(fctr)*fftshift(res(:,:,:,n));
    end
end

if ( numel(s) == 5 )
    for s=1:size(x,5)        
        for n=1:size(x,4)        
            res(:,:,:,n,s) = ifftshift(x(:,:,:,n,s));
            res(:,:,:,n,s) = ifftn(res(:,:,:,n,s));
            res(:,:,:,n,s) = sqrt(fctr)*fftshift(res(:,:,:,n,s));
        end
    end
end
