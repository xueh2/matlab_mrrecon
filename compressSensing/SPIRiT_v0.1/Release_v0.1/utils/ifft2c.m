function res = ifft2c(x)

s = size(x);
fctr = size(x,1)*size(x,2);

if ( numel(s) == 2 ) % Nfe Npe
    res = sqrt(fctr)*fftshift(ifft2(ifftshift(x)));
end

if ( numel(s) == 3 ) % Nfe Npe Frame   
    for n=1:size(x,3)
        res(:,:,n) = sqrt(fctr)*fftshift(ifft2(ifftshift(x(:,:,n))));
    end
end

if ( numel(s) == 4 ) % Nfe Npe Cha Frame
    for n=1:size(x,4)
        for c=1:size(x,3)
            res(:,:,c,n) = sqrt(fctr)*fftshift(ifft2(ifftshift(x(:,:,c,n))));
        end
    end
end

if ( numel(s) > 4 )
    for u=1:size(x,9)
        for r=1:size(x,8)
            for q=1:size(x,7)
                for p=1:size(x,6)
                    for t=1:size(x,5)
                        for s=1:size(x,4)
                            for n=1:size(x,3)
                                res(:,:,n,s,t,p,q,r,u) = sqrt(fctr)*fftshift(ifft2(ifftshift(x(:,:,n,s,t,p,q,r,u))));
                            end
                        end
                    end
                end
            end
        end
    end
end