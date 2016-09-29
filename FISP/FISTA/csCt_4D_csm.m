%
% un compressed-sensing function
%

function vICSimage = csCt(vCSimage)

global sizeImage;
% it contains five elements  [1d,2d,3d,c,t]

global samplingLocations;
% it contains four elements [1d,2d,3d,t]

global coilSensitivity;
% it contains four elements  [1d,2d,3d,c]
%          or five  elements [1d,2d,3d,c,t]

% vCSimage is of size (t x nnz(samplingLocations) ) x c

vICSimage=zeros(prod(sizeImage(1:3)),sizeImage(5));

switch numel(size(coilSensitivity))
    case 4
        index_e=0;
        for i=1:sizeImage(5)
            sLocations=samplingLocations(:,:,:,i);
            index_b=index_e+1;
            index_e=index_e+nnz(sLocations);
            
            for j=1:sizeImage(4)
                A=zeros(sizeImage(1:3));
                A(logical(sLocations))=vCSimage(index_b:index_e,j);
                
                ICSimage = conj(coilSensitivity(:, :, :,j)) .* ifftn(A) ...
                    * sqrt(prod(sizeImage(1:3)));
                
                vICSimage(:,i)=vICSimage(:,i)+ICSimage(:);
            end            
        end
        
    case 5
        index_e=0;
        for i=1:sizeImage(5)
            sLocations=samplingLocations(:,:,:,i);
            index_b=index_e+1;
            index_e=index_e+nnz(sLocations);
            
            for j=1:sizeImage(4)
                A=zeros(sizeImage(1:3));
                A(logical(sLocations))=vCSimage(index_b:index_e,j);
                
                ICSimage = conj(coilSensitivity(:, :, :,j,i)) .* ifftn(A) ...
                    * sqrt(prod(sizeImage(1:3)));
                
                vICSimage(:,i)=vICSimage(:,i)+ICSimage(:);
            end            
        end
        
    otherwise
        error();
end