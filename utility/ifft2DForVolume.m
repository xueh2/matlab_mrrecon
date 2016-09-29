
function Img = ifft2DForVolume(K_space, centreFFT)
% perform 2D ifft for volume

s = size(K_space);
if length(s) == 2, s(3) = 1; end
if (length(s) ~=3) & (length(s) ~=4)
    'Error! The input data must be 3D/4D k-space data!'
    return,
end

r = sqrt(size(K_space,1)*size(K_space,2));

if ( length(s) == 3 )
    Img = zeros( s(1), s(2), s(3) );

    if nargin == 1
        for z_index=1:s(3)
            Img(:,:,z_index) = r*ifftshift(ifft2(fftshift(K_space(:, :, z_index))));
        end
    elseif nargin == 2
        for z_index=1:s(3)
            Img(:,:,z_index) = r*ifft2(fftshift(K_space(:, :, z_index)));
        end
    end
end

if ( length(s) == 4 )
    Img = zeros( s(1), s(2), s(3), s(4) );

    if nargin == 1
        for f=1:s(4)
            for c=1:s(3)
                Img(:,:,c,f) = r*ifftshift(ifft2(fftshift(K_space(:, :, c, f))));
            end
        end
    elseif nargin == 2
        for f=1:s(4)
            for c=1:s(3)
                Img(:,:,c,f) = r*ifft2(fftshift(K_space(:, :, c, f)));
            end
        end
    end
end
