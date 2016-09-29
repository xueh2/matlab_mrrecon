function yy = soft_csm_4D_gs(x,T)

yy=zeros(size(x));

% applying L1 thresholding to each entry
y=max( abs(x) - T(1), 0 );
y=y./(y+T(1)).* x;

% applying L2 thresholding to each row
a=sqrt(sum( abs(y).^2, 2)); % compute the L2 norm of each row
b=max( a - T(2), 0); 
yy= repmat( a./(b + T(2)), 1, size(x,2)).*y;
