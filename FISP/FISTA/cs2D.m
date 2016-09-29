%
% compressed-sensing function
%

% for 2D cartesian sampling without CSM

function y = cs2D(vImage)

global sizeImage;
% it contains three elements [1d,2d,c] 

global samplingLocations;
% it is a 3D tensor of size [1d,2d, c]

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
    
       
    case 3
        
        imageRecon = reshape(vImage, sizeImage);
        
        y=[];
        for i=1:sizeImage(3)
            imageReconCS = fft2(imageRecon(:,:,i)) / sqrt(prod(sizeImage(1:2)));
            
            sLocations=samplingLocations(:,:,i);
            
            vCSimage = sLocations.* imageReconCS;
            
            yy = vCSimage(sLocations(:) == 1);
            y=[y; yy];
        end
            
    otherwise
           error('\n error in sizeImage');            
end
