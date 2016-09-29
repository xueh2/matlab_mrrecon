function y = softThresholding(x,mu)
%soft tresholding operator
    y=abs(x)-mu;
    y=(y+abs(y))/2;
    y=y.*sign(x);
end
