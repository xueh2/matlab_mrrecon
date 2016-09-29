function y = soft_temporal_2d(x,T, scale)

% x is temporal 3D

n=size(x,2);

y(:,1:(n/2)) = max(abs( x(:,1:(n/2)) ) - T, 0);
y = y(:,1:(n/2))./(y(:,1:(n/2))+T) .* x(:,1:(n/2));

y(:,(n/2+1):n) = max(abs( x(:,(n/2+1):n ) ) - scale*T, 0);
y(:,(n/2+1):n) = y(:,(n/2+1):n )./(y(:,(n/2+1):n )+scale*T) .* x(:,(n/2+1):n);

y=y(:);