
function kspace = fft2DForVolume(Im, centreFFT)
% perform 2D ifft for volume

s = size(Im);
if length(s) == 2, s(3) = 1; end
if (length(s) ~=3) & (length(s) ~=4)
    'Error! The input data must be 3D/4D complex image data!'
    return,
end

r = 1/sqrt(size(Im,1)*size(Im,2));

if ( length(s) == 3 )
    kspace = zeros( s(1), s(2), s(3) );

    if nargin == 1
        for z_index=1:s(3)
            kspace(:,:,z_index) = r*ifftshift(fft2(fftshift(Im(:, :, z_index))));
        end
    elseif nargin == 2
        for z_index=1:s(3)
            kspace(:,:,z_index) = r*fft2(fftshift(Im(:, :, z_index)));
        end
    end
end

if ( length(s) == 4 )
    kspace = zeros( s(1), s(2), s(3), s(4) );

    if nargin == 1
        for f=1:s(4)
            for c=1:s(3)
                kspace(:,:,c,f) = r*ifftshift(fft2(fftshift(Im(:, :, c, f))));
            end
        end
    elseif nargin == 2
        for f=1:s(4)
            for c=1:s(3)
                kspace(:,:,c,f) = r*fft2(fftshift(Im(:, :, c, f)));
            end
        end
    end
end
