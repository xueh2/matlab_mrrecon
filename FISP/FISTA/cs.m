%
% compressed-sensing function
%

function vCSimage = cs(vImage)

global sizeImage;
global samplingLocations;
global coilSensitivity;

imageRecon = reshape(vImage, sizeImage(1:2));
imageReconCS = zeros(sizeImage);
  
switch numel(sizeImage)
  case 2
    imageReconCS = fft2(imageRecon) / sqrt(prod(sizeImage(1:2)));
  case 3
    for slice2Dnum = 1:sizeImage(3)
      imageReconCS(:, :, slice2Dnum) = ...
        fft2(coilSensitivity(:, :, slice2Dnum) .* imageRecon) ...
        / sqrt(prod(sizeImage(1:2)));
    end
  otherwise
    error ();
end

vCSimage = samplingLocations .* imageReconCS;
vCSimage = vCSimage(samplingLocations(:) == 1);
