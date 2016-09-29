
function Img = SoS_Image_TemporalArray(Im)
% Sum of Square image reconstruction

s = size(Im);
ImgArray = zeros( s(1), s(2), s(4) );
for f=1:s(4)
    aIm = SoS_Image(Im(:,:,:,f));
    ImgArray(:,:,f) = aIm;
end
Img = ImgArray;
