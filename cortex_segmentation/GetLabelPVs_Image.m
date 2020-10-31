
function LabelPVs_Image = GetLabelPVs_Image(mix)
% transform posterior as image file

xsize = mix.header.xsize;
ysize = mix.header.ysize;
zsize = mix.header.zsize;

LabelPVs_Image = zeros([ysize xsize zsize], 'uint32');
if ( isempty(mix.LabelPVs) == 0 )
    LabelPVs_Image(mix.LabelPVs(:)) = 1;
end
return;