function y = soft_2D(x,w,T)

% x is a 2D tensor, and the penalty is |x|_{2,1}
% w is a 1D vector
% T is a scalar

y = max( sqrt( sum(abs(x).^2,2) ) - w.*T, 0);
y = x.* repmat(y./(y+ w*T), 1, size(x,2));
