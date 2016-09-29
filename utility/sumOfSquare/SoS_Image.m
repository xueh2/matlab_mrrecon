
% Sum of Square image reconstruction
function Img = SoS_Image(im)

s = size(im);

if ( numel(s) == 2 )
    Img = abs(im);
end

if ( numel(s) == 3 )
    Img = zeros(s(1), s(2));
    for z_index=1:s(3)
        ima = im(:, :, z_index);
        Img = Img + conj(ima).*ima;
    end
    Img = sqrt(Img);
end

if ( numel(s) == 4 )
    Img = zeros(s(1), s(2), s(3));
    for z_index=1:s(4)
        ima = im(:, :, :, z_index);
        Img = Img + conj(ima).*ima;
    end
    Img = sqrt(Img);
end
