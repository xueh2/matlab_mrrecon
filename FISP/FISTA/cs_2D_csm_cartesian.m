%
% compressed-sensing function
%
% vCSimage=A * x

% for Cartesian Sampling

function vCSimage = cs_2D_csm(vImage)

global sizeImage;
% it contains five elements  [1d,2d,c,t]

global samplingLocations;
% it contains four elements [1d,2d,t]

global coilSensitivity;
% it contains four elements  [1d,2d,c]
%          or five  elements [1d,2d,c,t]

% vImage is of size: (1d x 2d) x t
% vCSimage is of size: sampled (1d x 2d x t) x # channel

imageRecon = reshape(vImage, [sizeImage(1:2), sizeImage(4)]);

vCSimage=zeros(nnz(samplingLocations), sizeImage(3));

switch numel(size(coilSensitivity))
    case 3
        index_e=0;
        for i=1:sizeImage(4)
            sLocations=samplingLocations(:,:,i);
            index_b=index_e+1;
            index_e=index_e+nnz(sLocations);
            
            for j=1:sizeImage(3)
                imageReconCS = fft2(coilSensitivity(:,:,j).*imageRecon(:,:,i)) / sqrt(prod(sizeImage(1:2)));
                
                vCSimage(index_b:index_e,j)=imageReconCS(logical(sLocations(:)));
            end
        end
        
    otherwise
        error();
end







