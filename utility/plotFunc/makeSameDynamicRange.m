
function adjIm = makeSameDynamicRange(base, im)
% make the dynamic range of im similar to the base by aligning the median value

i = median(base(:))
o = median(im(:))

if ( i ~= 0 )
    r = o/i;
else
    r = 1;
end

adjIm = im;
adjIm = adjIm / r;
