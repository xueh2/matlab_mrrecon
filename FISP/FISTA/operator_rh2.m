function Y=operator_rh2(x, w,th, W, WT)

global sizeImage 
global sizeReconImage;
global analysisSize;

switch numel(sizeImage)
    case 3 % for 2D data
        for i=1:size(x,2)
            X(:,i)=soft( W( x(:,i) ), th(i));
        end
        
        for i=1:size(x,2)
            Y(:,i)=WT( X(:,i) );
        end
        
    case 5 % for 4D MRA data            % [112, 84, 30, 4, 13];
        % x(:,1:4,1:13)
        for i=1:sizeImage(5)
            V=[];
            index=1:prod(sizeImage(1:2));
            for k=1:sizeImage(3)
                X=[];
                for j=1:sizeImage(4)
                    X(:,j)=soft( W( x(index,j,i) ), th(i) );
                end
                
                U=[];
                for j=1:sizeImage(4)
                    U(:,j)=WT( X(:,j) );
                end
                V=[V; U];
                
                index=index+prod(sizeImage(1:2));
            end
            Y(:,:,i)=V;
        end
        
end