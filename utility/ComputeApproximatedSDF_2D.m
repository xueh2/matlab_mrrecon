function segSDF = ComputeApproximatedSDF_2D(seg2D)
% for a 2D+T input segmenation array, compute corresponding signed distance field

N = size(seg2D, 3);

segSDF = zeros(size(seg2D));

for n=1:N
    s2D = seg2D(:,:,n);
    s2D = double(s2D);
    
    ind = find(s2D(:)>0);
    if (~isempty(ind))
        sdf2D = CreateApproximated_SDF_2D(s2D);       
        segSDF(:,:,n) = sdf2D;
    end
end
