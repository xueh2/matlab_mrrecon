
function Img = SoS_TemporalArray(K_space, c)
% Sum of Square image reconstruction
% K_psace: [RO E1 CHA N]

s = size(K_space);
ImgArray = zeros( s(1), s(2), s(4) );
for f=1:s(4)
    Img = zeros(s(1), s(2));
    if nargin == 1   
        for z_index=1:s(3)
            % ima = ifftshift(ifft2(ifftshift(K_space(:, :, z_index, f))));
            ima = ifft2c(K_space(:, :, z_index, f));
            Img = Img + conj(ima).*ima;
        end
    elseif nargin == 2
        for z_index=1:s(3)
            ima = ifft2(ifftshift(K_space(:, :, z_index, f)));
            % ima = ifft2c(K_space(:, :, z_index, f));
            Img = Img + conj(ima) .* ima;
        end
    end
    Img = sqrt(Img);
    ImgArray(:,:,f) = Img;
end
Img = ImgArray;
