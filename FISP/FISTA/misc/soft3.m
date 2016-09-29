function y = soft3(x,w,T)

% x is a matrix, and the penalty is |x|_{2,1}

global sizeImage;

switch numel(sizeImage)
       
    case 5
        for i=1: size(x,3)
            yy = max( sqrt( sum(abs(x(:,:,i)).^2,2) ) - w(:,i).*T(i), 0);
            yy = x(:,:,i).* repmat(yy./(yy+ w(:,i).*T(i)), 1, size(x,2));
            
            y(:,:,i)=yy;
        end        
        
    otherwise
        error('check x');
end