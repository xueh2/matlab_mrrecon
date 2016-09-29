function y = soft_4D(x,w,T,pFlag)

% x is a 3D tensor, and the penalty is |x|_{2,1}
% w is a 2D matrix
% T is a 1D vector

% global sizeImage;
% 
% if numel(sizeImage)~=5
%     error('\n The input parameter sizeImage should have 5 elements');
% end

if pFlag==1
    for i=1: size(x,3)
        yy = max( sqrt( sum(abs(x(:,:,i)).^2,2) ) - w(:,i).*T(i), 0);
        yy = x(:,:,i).* repmat(yy./(yy+ w(:,i).*T(i)), 1, size(x,2));
        
        y(:,:,i)=yy;
    end
    
end

if pFlag==4
    for i=1: size(x,3)
        yy = max( abs(x(:,:,i))-T(i), 0);
        yy = x(:,:,i).* (yy./(yy+ T(i)));
        
        y(:,:,i)=yy;
    end
end
        