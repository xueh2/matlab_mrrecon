function w = hanning1d(n);
%hanning1d   Symmetric Hanning window. 


if ~rem(n,2); % n even
%     m=n/2;
%     for i=1:m
%         w(i) = .5*(1 - cos(2*pi*i/(n+1))); 
%     end
%     for i=m+1:n
%         w(i)=w(n+1-i);
%     end
%     
%     
    for i=0:n-1
        w(i+1) = .5*(1 - cos(2*pi*(i+1)/(n+1))); 
    end
    
    
else % n odd
    m = (n+1)/2;
%     for i=1:m
%         w(i) = .5*(1 - cos(2*pi*i/(n+1))); 
%     end
%     for i=m+1:n
%         w(i)=w(n+1-i);
%     end
    
    for i=0:m-1
        w(i+1) = .5*(1 - cos(2*pi*(i+1)/(n+1))); 
    end
    for i=m:n-1
        w(i+1)=w(n-i);
    end
    
    
    
end

   