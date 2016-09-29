function [volumeResampled, headerResampled] = resampleVolume(srcVolume, headerSrc, dstVolume, headerDst)
% resample srcVolume on the grid defined by dstVolume
% the BSpline coefficient is computed and BSpline interpolation is used

% compute BSpline coefficient
coeff = Matlab_ComputeBSplineCoefficient(double(srcVolume), headerSrc, 'Row-wise');
[volumeResampled, headerResampled] = Matlab_EvaluateBSplineResampledSlices(double(srcVolume), coeff, headerSrc, double(dstVolume), headerDst, 'BSpline', 'Row-wise');

% % % for every pixel point in the dstVolume
% volumeResampled = zeros(size(dstVolume));
% headerResampled = headerDst;
% 
% for z=1:headerDst.sizeZ
%     for y=1:headerDst.sizeY
%         for x=1:headerDst.sizeX
%             % image to world
%             [wx, wy, wz] = Image2WorldMrFtk(headerDst, x-1, y-1, z-1);            
%             % world to src image
%             [ix, iy, iz] = World2ImageMrFtk(headerSrc, wx, wy, wz);            
%             % evaluate the pixel value
%             volumeResampled(y, x, z) = Matlab_EvaluateBSplineInterpolation(coeff, ix, iy, iz, 'Row-wise');
%         end
%     end
% end
            