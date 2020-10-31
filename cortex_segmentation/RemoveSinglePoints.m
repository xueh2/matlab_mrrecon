
function suspected_volume = RemoveSinglePoints(suspected_volume, volumeThres)
% remove single points (outlier ... ?)

[l,num] = bwlabeln(suspected_volume, 6);
volumes=zeros(num,1);

total = numel(suspected_volume);
for k=1:total
    if (l(k)>0)
        volumes(l(k))=volumes(l(k))+1;
    end
end

for k = 1:total
    if (l(k)>0)
        if ( (volumes(l(k)) <= volumeThres ) )
            suspected_volume(k) = 0;
        end
    end
end
return