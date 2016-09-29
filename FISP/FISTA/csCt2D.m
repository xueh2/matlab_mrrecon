%
% un compressed-sensing function
%

% for 2D cartesian sampling without CSM

function y = csCt2D(vCSimage)

global sizeImage;
% it contains three elements [1d,2d,c] 

global samplingLocations;
% it is a 3D tensor of size [1d,2d,c]

%global coilSensitivity;

% temp = zeros(prod(sizeImage), 1);
% temp(samplingLocations(:) ~= 0) = vCSimage;
% vCSimage = temp;
%   
% CSimage = reshape(vCSimage, sizeImage);
% ICSimage = zeros(sizeImage(1:2));
%
% switch numel(sizeImage)
%   case 2
%     ICSimage = ifft2(samplingLocations .* CSimage) * sqrt(prod(sizeImage(1:2)));
%   case 3
%     for slice2Dnum = 1:sizeImage(3)
%       ICSimage = ICSimage + ...
%         conj(coilSensitivity(:, :, slice2Dnum)) .* ifft2(CSimage(:, :, slice2Dnum)) ...
%         * sqrt(prod(sizeImage(1:2)));
%     end
%   otherwise
%     error();
% end
% 
% vICSimage = ICSimage(:);


switch numel(sizeImage)
            
    case 3
        
        index_e=0;
        for i=1:sizeImage(3)
            index_b=index_e+1;
            index_e=index_e+ nnz(samplingLocations(:,:,i));
            
            temp=zeros(prod(sizeImage(1:2)), 1);
            
            sLocations=samplingLocations(:,:,i);
            
            temp(sLocations(:) ~= 0) = vCSimage(index_b:index_e,1);
            
            CSimage = reshape(temp, sizeImage(1:2));
            
            ICSimage = ifft2(CSimage) * sqrt(prod(sizeImage(1:2)));
            
            y(:,i)=ICSimage(:);
        end
        
    otherwise
        error('\n error in sizeImage');
end
