function [endo, epi, rv, rvi] = Seg2contours(Seg, new_len, n_components)

    if (nargin<2)
        new_len = 200;
    end

    if (nargin<3)
        n_components = 50;
    end
    
    SLC = numel(Seg);
    
    endo = cell(SLC, 1);
    epi = cell(SLC, 1);
    rv = cell(SLC, 1);
    rvi = zeros(SLC, 2);
    
    for s=1:SLC
        endo{s} = mask2contour(Seg(s).endo_resized_mask_norm, 1, new_len, n_components);
        epi{s} = mask2contour(Seg(s).epi_resized_mask_norm, 1, new_len, n_components);
        rv{s} = mask2contour(Seg(s).rv_resized_mask_norm, 1, new_len, n_components);
        
        [I, J] = find(Seg(s).rvi_resized_mask_norm>0);        
        rvi(s, 1) = mean(I);
        rvi(s, 2) = mean(J);
    end
end