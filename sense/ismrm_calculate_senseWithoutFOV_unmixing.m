function [unmix, gmap] = ismrm_calculate_senseWithoutFOV_unmixing(acc_factor, csm, snr_map, alpha)
%
% function [unmix, gmap] = ismrm_calculate_senseWithoutFOV_unmixing(acc_factor, csm)
%
% Calculates the unmixing coefficients for a 2D image
%
% INPUT:
%       acc_factor  scalar       : Acceleration factor, e.g. 2
%       csm         [x, y, coil, 2] : Coil sensitivity map
%
% OUTPUT:
%       unmix       [x, y, coil, 2] : Image unmixing coefficients for a single x location 
%       gmap        [x, y]       : Noise enhancement map 
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip Beatty (philip.beatty@sri.utoronto.ca)
%

unmix = zeros(size(csm));

for x=1:size(csm,1), 
    unmix(x,:,:,:) = ismrm_calculate_senseWithoutFOV_unmixing_1d(acc_factor, squeeze(csm(x,:,:,:)), squeeze(snr_map(x, :)), alpha); 
end

if (nargout > 1),
   gmap = sqrt(sum(abs(unmix(:,:,:,1)).^2,3)).*sqrt(sum(abs(csm(:,:,:,1)).^2,3));
end
