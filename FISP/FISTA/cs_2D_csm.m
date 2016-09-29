%
% compressed-sensing function
%
% vCSimage=A * x


function vCSimage = cs_2D_csm(vImage)

global sizeImage;
% it contains five elements  [1d,2d,c,t]

global coilSensitivity;
% it contains four elements  [1d,2d,c]
%          or five  elements [1d,2d,c,t]

global kcoor;
% the k-space coordinate
% a 3D matrix [# samples, 2, t]
% where # samples denote the number of sampled locations in k-space
% kcoor(:,1)-- x-axis
% kcoor(:,2)-- y-axis
% t denotes the number of temporal phases

% vImage is of size (1d x 2d) x t

imageRecon = reshape(vImage, [sizeImage(1:2), sizeImage(4)]);

vCSimage=zeros( size(kcoor,1), sizeImage(3));

switch numel(size(coilSensitivity))
    case 3
        index_e=0;
        for i=1:sizeImage(4)
            index_b=index_e+1;
            index_e=index_e+ size(kcoor,1);
            
            ik=kcoor(:,1,i)/(2*pi)+1i*kcoor(:,2,i)/(2*pi); w = 1; phase = 1; shift = [0,0];
            nuFt=NUFFT(ik,w,shift,sizeImage(1:2));

            for j=1:sizeImage(3)
                imageReconCS = nuFt * (coilSensitivity(:,:,j).*imageRecon(:,:,i)); 
                
                vCSimage(index_b:index_e,j)=imageReconCS;
            end
        end
        
    case 4
        
        index_e=0;
        for i=1:sizeImage(4)
            index_b=index_e+1;
            index_e=index_e+ size(kcoor,1);
            
            ik=kcoor(:,1,i)/(2*pi)+1i*kcoor(:,2,i)/(2*pi); w = 1; phase = 1; shift = [0,0];
            nuFt=NUFFT(ik,w,shift,sizeImage(1:2));
            
            for j=1:sizeImage(3)
                imageReconCS = nuFt * (coilSensitivity(:,:,j,i).*imageRecon(:,:,i)) / sqrt(prod(sizeImage(1:2)));
                
                vCSimage(index_b:index_e,j)=imageReconCS;
            end            
        end
        
    otherwise
        error();
end

