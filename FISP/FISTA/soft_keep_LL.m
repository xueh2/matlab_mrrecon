function y = soft_keep_LL(x,T)

y = max(abs(x) - T, 0);
y = y./(y+T) .* x;

N = numel(x);
y(1:N/2) = x(1:N/2);