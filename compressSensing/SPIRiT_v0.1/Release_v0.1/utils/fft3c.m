function res = fft3c(x)
fctr = size(x,1)*size(x,2)*size(x,3);

s = size(x);

res = x;

if ( length(s) == 3 )   
    res = ifftshift(x);   
    res = fftn(res);   
    res = 1/sqrt(fctr)*fftshift(res);
end

if ( length(s) == 4 )
    for s=1:size(x,4)
        res(:,:,:,s) = ifftshift(x(:,:,:,s));
        res(:,:,:,s) = fftn(res(:,:,:,s));
        res(:,:,:,s) = 1/sqrt(fctr)*fftshift(res(:,:,:,s));
    end
end

if ( length(s) == 5 )
    for s5=1:size(x,5)
        for s=1:size(x,4)
            res(:,:,:,s,s5) = ifftshift(x(:,:,:,s,s5));
            res(:,:,:,s,s5) = fftn(res(:,:,:,s,s5));
            res(:,:,:,s,s5) = 1/sqrt(fctr)*fftshift(res(:,:,:,s,s5));
        end
    end
end


