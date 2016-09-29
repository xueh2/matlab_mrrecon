%
% un compressed-sensing function
%

% for Cartesian Sampling

function vICSimage = csCt_2D_csm(vCSimage)

global sizeImage;
% it contains five elements  [1d,2d,c,t]

global coilSensitivity;
% it contains four elements  [1d,2d,c]
%          or five  elements [1d,2d,c,t]

global samplingLocations;
% it contains four elements [1d,2d,t]

% vCSimage is of size: sampled (1d x 2d x t) x c
% vICSimage is of size: (1d x 2d) x t

vICSimage=zeros(prod(sizeImage(1:2)),sizeImage(4));

switch numel(size(coilSensitivity))
    case 3
        index_e=0;
        for i=1:sizeImage(4) % t
            sLocations=samplingLocations(:,:,i);
            index_b=index_e+1;
            index_e=index_e+nnz(sLocations);
            
            for j=1:sizeImage(3) % c
                A=zeros(sizeImage(1:2));
                A(logical(sLocations))=vCSimage(index_b:index_e,j);
                
                ICSimage = conj(coilSensitivity(:, :, j)) .* ifft2(A) ...
                    * sqrt(prod(sizeImage(1:2)));
                
                vICSimage(:,i)=vICSimage(:,i)+ICSimage(:);
            end            
        end

        
    otherwise
        error();
end

