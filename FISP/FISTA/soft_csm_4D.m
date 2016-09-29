function yy = soft_csm_4D(x,T)

yy=zeros(size(x));
for i=1:size(x,2)
    y = max(abs(x(:,i)) - T(i), 0);
    y = y./(y+T(i)) .* x(:,i);
    yy(:,i)=y;
end
