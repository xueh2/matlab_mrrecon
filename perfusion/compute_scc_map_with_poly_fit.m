function [scc_map_final, scc_map, polyfit_output] = compute_scc_map_with_poly_fit(pd, mask, bandwidth)
% [scc_map_final, scc_mask, polyfit_output] = compute_scc_map_with_poly_fit(pd, bandwidth)

    degree=2;

    threshold_mask = mask;
    [Ny,Nx]=size(pd);
    [x,y]=meshgrid(1:Nx,1:Ny);

    p=polyfit2dw(x,y,pd,mask,degree,degree);
    polyfit_output=polyval2d(p,x,y);

    [scc_map, soheilmask]=Matlab_gt_surface_coil_correction(pd, pd,[],0,11,3.0,1);

    % combine 
    [Y, X] = find(mask>0);
    ns = createns([X, Y],'nsmethod','kdtree');
    
    [Y2, X2] = find(mask==0);
    [idx, dist] = knnsearch(ns, [X2 Y2], 'k', 1);

    dist_map = zeros(size(mask));

    for k=1:size(Y2)
        dist_map(Y2(k), X2(k)) = dist(k);
    end

    weights_map = ones(size(mask));

    ind = find(mask>0);
    mean_pd = mean(pd(ind))

    RO = size(mask, 1);
    E1 = size(mask, 2);
    for e1=1:E1
        for ro=1:RO
            if(dist_map(ro, e1)==0)
                weights_map(ro,e1) = 0;
            elseif(dist_map(ro, e1)<=bandwidth)
    %             if(scc(ro,e1)>output(ro,e1))
    %                 weights_map(ro,e1) = 0.5;
    %             else
                    weights_map(ro,e1) = 0.0;
    %             end
                    if(pd(ro,e1)>0.35*mean_pd)
                        weights_map(ro,e1) = 1.0;
                    end
%                     if(scc_map(ro,e1)>polyfit_output(ro,e1))
%                         weights_map(ro,e1) = 1.0;
%                     end
            else
                t = 5*(dist_map(ro, e1)-bandwidth)/bandwidth;
                weights_map(ro,e1) = (1+tanh(t))/2;
                %weights_map(ro,e1) = 1;
                if(pd(ro,e1)>0.1*mean_pd)
                    weights_map(ro,e1) = 1.0;
                end
            end

        end
    end

    % figure; imagescn(weights_map);
    
    scc_map_final = weights_map.*scc_map + (1.0-weights_map).*polyfit_output;

    ind = find(isnan(scc_map_final) | scc_map_final<1);
    if(numel(ind)>0)
        scc_map_final(ind) = scc_map(ind);
    end
    
    mean_scc_level = mean(scc_map_final(ind))
    
    for e1=1:E1
        for ro=1:RO
            if(dist_map(ro, e1)>0)
                if(scc_map(ro,e1)>scc_map_final(ro,e1))
                    scc_map_final(ro,e1) = scc_map(ro,e1);                    
                end
                if(scc_map_final(ro,e1)<0.5*mean_scc_level)
                    scc_map_final(ro,e1) = mean_scc_level;                    
                end
            end
        end
    end
    
    scc_map_final = medfilt2(scc_map_final, [round(bandwidth/2) round(bandwidth/2)]);
end

