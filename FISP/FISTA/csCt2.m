%
% un compressed-sensing function
%

function y = csCt2(vCSimage)

global sizeImage;
% it contains four elements [1d,2d,3d,c] 
%         or five elements  [1d,2d,3d,c,t]

global samplingLocations;
% it is a 3D tensor of size [1d,2d,3d]
%      or 4D tensor of size [1d,2d,3d,t]

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
    
    case 4
        
        for i=1:sizeImage(4)
            temp=zeros(prod(sizeImage(1:3)), 1);
            temp(samplingLocations(:) ~= 0) = vCSimage(:,i);
            
            CSimage = reshape(temp, sizeImage(1:3));
            
            ICSimage = ifftn(CSimage) * sqrt(prod(sizeImage(1:3)));
            
            y(:,i)=ICSimage(:);
        end
        
        
    case 5
        
        for j=1:sizeImage(5)
            for i=1:sizeImage(4)
                temp=zeros(prod(sizeImage(1:3)), 1);
                
                if (numel(size(samplingLocations)) ==3)
                    temp(samplingLocations(:) ~= 0) = vCSimage(:,i,j);
                else
                    if (numel(size(samplingLocations)) ==4)
                        sLocations=samplingLocations(:,:,:,j);
                        temp(sLocations(:) ~= 0) = vCSimage(:,i,j);
                    end
                end
                
                CSimage = reshape(temp, sizeImage(1:3));
                
                ICSimage = ifftn(CSimage) * sqrt(prod(sizeImage(1:3)));
                
                y(:,i,j)=ICSimage(:);
            end
        end
        
    otherwise
        error('\n error in sizeImage');
end
