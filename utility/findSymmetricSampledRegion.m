function range = findSymmetricSampledRegion(rangeSampled, centerInd)
% find the symmetric sampling range

halfSizeStart = centerInd - rangeSampled(1);
halfSizeEnd =  rangeSampled(2) - centerInd;

if ( halfSizeStart > halfSizeEnd )
    range(1) = centerInd - halfSizeEnd;
    range(2) = centerInd + halfSizeEnd;
else
    range(1) = centerInd - halfSizeStart;
    range(2) = centerInd + halfSizeStart;
end