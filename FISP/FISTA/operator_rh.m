function Y=operator_rh(x, w,th, W, WT)

global sizeImage 
global sizeReconImage;
global analysisSize;

switch numel(sizeImage)
    
    case 3 % for 2D data
        
        for i=1:size(x,2)
            V(:,i)=W( x(:,i) );
        end
        X=soft_2D(V, 1, th);
        
        for i=1:size(x,2)
            Y(:,i)=WT( X(:,i) );
        end
        
    case 5 % for 4D MRA data            % [112, 84, 30, 4, 13];
           % x(:,1:4,1:13)
        for i=1:sizeImage(5)
            V=[];
            index=1:prod(sizeImage(1:2));
            for k=1:sizeImage(3)
                U=[];
                for j=1:sizeImage(4)
                    U(:,j)=W( x(index,j,i) );
                end
                X=soft_2D(U, 1, th(i));
                % X is the result for a 2D image
                
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