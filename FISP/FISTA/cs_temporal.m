%
% compressed-sensing function
%

% for 2D cartesian sampling with CSM

function vCSimage = cs(vImage)

global sizeImage;
global samplingLocations;
global coilSensitivity;

imageRecon = reshape(vImage, sizeImage(1:2));
imageReconCS = zeros(sizeImage);
  
switch numel(sizeImage)
  case 3
    for slice2Dnum = 1:sizeImage(3)
      imageReconCS(:, :, slice2Dnum) = ...
        fft (coilSensitivity(:, :, slice2Dnum) .* imageRecon, [], 1) ...
        / sqrt(prod(sizeImage(1)));
    end
  otherwise
    error ();
end

vCSimage = samplingLocations .* imageReconCS;
vCSimage = vCSimage(samplingLocations(:) == 1);
