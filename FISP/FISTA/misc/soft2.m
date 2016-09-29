function y = soft2(x,T)

% x is a matrix, and the penalty is |x|_{2,1}

global sizeImage;

switch numel(sizeImage)
    case 4
        
        y = max( sqrt( sum(abs(x).^2,2) ) - T, 0);
        y = x.* repmat(y./(y+T), 1, size(x,2));
        
    case 5
        for i=1: size(x,3)
            yy = max( sqrt( sum(abs(x(:,:,i)).^2,2) ) - T(i), 0);
            yy = x(:,:,i).* repmat(yy./(yy+T(i)), 1, size(x,2));
            
            y(:,:,i)=yy;
        end        
        
    otherwise
        error('check x');
end