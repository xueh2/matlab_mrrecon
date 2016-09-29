
% coil reconstruction using sensitivity

function Img = SensitivityCoilCombination(Im, sensitivityMap)

s = size(Im);
sSen = size(sensitivityMap);
if length(s) == 3
    Img = zeros( s(1), s(2) );
    for c=1:s(3)
        Img = Img + conj(sensitivityMap(:,:,c)).*Im(:,:,c);
    end
end

if length(s) == 4
    Img = zeros( s(1), s(2), s(4) );
    for f=1:s(4)

        if ( length(sSen) == 4 )
            sen = sensitivityMap(:,:,:,f);
            for c=1:s(3)
                Img(:,:,f) = Img(:,:,f) + conj(sen(:,:,c)).*Im(:,:,c,f);
            end
        end

        if ( length(sSen) == 3 )
            for c=1:s(3)
                Img(:,:,f) = Img(:,:,f) + conj(sensitivityMap(:,:,c)).*Im(:,:,c,f);
            end
        end
    end
end





