%
% compressed-sensing function
%
% vCSimage=A * x


function vCSimage = cs(vImage)

global sizeImage;
% it contains five elements  [1d,2d,3d,c,t]

global samplingLocations;
% it contains four elements [1d,2d,3d,t]

global coilSensitivity;
% it contains four elements  [1d,2d,3d,c]
%          or five  elements [1d,2d,3d,c,t]

% vImage is of size (1d x 2d x 3d) x t

imageRecon = reshape(vImage, [sizeImage(1:3), sizeImage(5)]);

vCSimage=zeros(nnz(samplingLocations), sizeImage(4));

switch numel(size(coilSensitivity))
    case 4
        index_e=0;
        for i=1:sizeImage(5)
            sLocations=samplingLocations(:,:,:,i);
            index_b=index_e+1;
            index_e=index_e+nnz(sLocations);
            
            for j=1:sizeImage(4)
                imageReconCS = fftn(coilSensitivity(:,:,:,j).*imageRecon(:,:,:,i)) / sqrt(prod(sizeImage(1:3)));
                
                vCSimage(index_b:index_e,j)=imageReconCS(logical(sLocations(:)));
            end
        end
        
    case 5
        
        index_e=0;
        for i=1:sizeImage(5)
            sLocations=samplingLocations(:,:,:,i);
            index_b=index_e+1;
            index_e=index_e+nnz(sLocations);
            
            for j=1:sizeImage(4)
                imageReconCS = fftn(coilSensitivity(:,:,:,j,i).*imageRecon(:,:,:,i)) / sqrt(prod(sizeImage(1:3)));
                
                vCSimage(index_b:index_e,j)=imageReconCS(logical(sLocations(:)));
            end
        end
        
    otherwise
        error();
end

