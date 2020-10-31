function [best_pair, best_diff] = find_A_B_pair(A, B, pt_ids_A, pt_ids_B, GLS_4CH_A, GLS_4CH_B)
    ind_A = find(pt_ids_A==A);
    ind_B = find(pt_ids_B==B);    
    
    best_pair = [];
    best_diff = 1e3;
    for a=1:numel(ind_A)
        va = GLS_4CH_A(ind_A(a));
        if(va<0)
            continue;
        end
        for b=1:numel(ind_B)
            vb = GLS_4CH_B(ind_B(b));
            if(vb<0)
                continue;
            end 
            if(abs(va-vb)<best_diff)
                best_diff = abs(va-vb);
                best_pair = [ind_A(a), ind_B(b)];
            end
        end
    end
end