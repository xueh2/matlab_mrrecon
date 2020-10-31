function rsq = compute_rsquared(x,y,p)
% Computes the R^2 value of a fit given two vectors with the estimate (from
% fit) and measured data points

% p = polyfit(x,y,1)  
% yfit =  p(1) * x + p(2); % equivalently: yfit = polyval(p,x);

p = polyfit(x, y,1);

yfit =  p(1) * x + p(2);
SSresid = sum((y - yfit).^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal;

end
