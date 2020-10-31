function [GLS, MAPSE, LV_Lengh_PT1, LV_Lengh_PT2, L] = compute_GLS(pts, ES_as_first_phase)

    if(~isempty(find(pts<2)))
        GLS = -1;
        MAPSE = -1;
        LV_Lengh_PT1 = -1;
        LV_Lengh_PT2 = -1;
        L = -1;
        return
    end

    GLS = -1;
    MAPSE = -1;
    
    ptc = 0.5 * (pts(1,:,:) + pts(2,:,:));
    ptc = squeeze(ptc);
    
    N = size(pts, 3);
    
    L = zeros(N,1);
    LV_Lengh_PT1 = zeros(N, 1);
    LV_Lengh_PT2 = zeros(N, 1);
    
    for j = 1:N        
        L(j) = norm([ptc(1,j), ptc(2,j)] - [pts(3,1,j), pts(3,2,j)]);
        LV_Lengh_PT1(j) = norm([pts(1,1,j), pts(1,2,j)] - [pts(3,1,j), pts(3,2,j)]);
        LV_Lengh_PT2(j) = norm([pts(2,1,j), pts(2,2,j)] - [pts(3,1,j), pts(3,2,j)]);
    end

    if(ES_as_first_phase)
        GLS = 100 * (L(1) - L) / L(1);
    else
        maxL = max(L(:));
        sorted_L = sort(L);
        %GLS = 100 * (maxL - L) / mean(sorted_L(N-5:N));
        % maxL = (0.5 * (L(1)+L(end)));
        GLS = 100 * (maxL - L) / (0.5 * (L(1)+L(end)));
    end
    
    MAPSE = zeros(N, 1);
    ptc_ED = (ptc(:,1) + ptc(:,end)) / 2;
    for j = 1:N        
        MAPSE(j) = norm([ptc(1,j), ptc(2,j)] - [ptc_ED(1,1), ptc_ED(2,1)]);
    end
end