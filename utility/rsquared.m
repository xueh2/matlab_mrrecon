function rsq = rsquared(y,yfit)

% function rsq = rsquared(y,yfit)

% 

% Computes the R^2 value between 2 vectors
 
 
SSresid = sum((y - yfit).^2);

SStotal = (length(y)-1) * var(y);

rsq = 1 - SSresid/SStotal;
 
end

 