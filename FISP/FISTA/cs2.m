%
% compressed-sensing function
%

function y = cs2(vImage)

global sizeImage;
% it contains four elements [1d,2d,3d,c] 
%         or five elements  [1d,2d,3d,c,t]

global samplingLocations;
% it is a 3D tensor of size [1d,2d,3d]
%      or 4D tensor of size [1d,2d,3d,t]

%global coilSensitivity;

  
% switch numel(sizeImage)
%   case 2
%     imageReconCS = fft2(imageRecon) / sqrt(prod(sizeImage(1:2)));
%   case 3
%     for slice2Dnum = 1:sizeImage(3)
%       imageReconCS(:, :, slice2Dnum) = ...
%         fft2(coilSensitivity(:, :, slice2Dnum) .* imageRecon) ...
%         / sqrt(prod(sizeImage(1:2)));
%     end
%   otherwise
%     error ();
% end

switch numel(sizeImage)
    
    case 4        
        imageRecon = reshape(vImage, sizeImage);
        
        for i=1:sizeImage(4)
            imageReconCS = fftn(imageRecon(:,:,:,i)) / sqrt(prod(sizeImage(1:3)));
            
            vCSimage = samplingLocations .* imageReconCS;
            y(:,i) = vCSimage(samplingLocations(:) == 1);
        end
        
    case 5
        
        imageRecon = reshape(vImage, sizeImage);
        
        for j=1:sizeImage(5)
            for i=1:sizeImage(4)
                imageReconCS = fftn(imageRecon(:,:,:,i,j)) / sqrt(prod(sizeImage(1:3)));
                
                if (numel(size(samplingLocations)) ==3)                  
                    vCSimage = samplingLocations .* imageReconCS;
                    y(:,i,j) = vCSimage(samplingLocations(:) == 1);
                else
                    if (numel(size(samplingLocations)) ==4)
                        
                        sLocations=samplingLocations(:,:,:,j);
                        
                        vCSimage = sLocations.* imageReconCS;
                        
                        y(:,i,j) = vCSimage(sLocations(:) == 1);
                    end
                end                
            end
        end
            
    otherwise
           error('\n error in sizeImage');            
end
