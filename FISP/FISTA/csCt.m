%
% un compressed-sensing function
%

function vICSimage = csCt(vCSimage)

global sizeImage;
global samplingLocations;
global coilSensitivity;

temp = zeros(prod(sizeImage), 1);
temp(samplingLocations(:) ~= 0) = vCSimage;
vCSimage = temp;
  
CSimage = reshape(vCSimage, sizeImage);
ICSimage = zeros(sizeImage(1:2));

switch numel(sizeImage)
  case 2
    ICSimage = ifft2(samplingLocations .* CSimage) * sqrt(prod(sizeImage(1:2)));
  case 3
    for slice2Dnum = 1:sizeImage(3)
      ICSimage = ICSimage + ...
        conj(coilSensitivity(:, :, slice2Dnum)) .* ifft2(CSimage(:, :, slice2Dnum)) ...
        * sqrt(prod(sizeImage(1:2)));
    end
  otherwise
    error();
end

vICSimage = ICSimage(:);
