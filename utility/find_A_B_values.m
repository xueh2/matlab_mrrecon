function [ra, rb] = find_A_B_values(A, B, pt_ids_A, pt_ids_B, GLS_4CH_A, GLS_4CH_B)
    ind_A = find(pt_ids_A==A);
    ind_B = find(pt_ids_B==B);    
    
    ra = [];
    rb = [];
    for a=1:numel(ind_A)
        va = GLS_4CH_A(ind_A(a));
        if(va>0 & va<30)
            ra = [ra va];
        end        
    end
    for b=1:numel(ind_B)
        vb = GLS_4CH_B(ind_B(b));
        if(vb>0 & vb<30)
            rb = [rb vb];
        end 
    end
    
    if(numel(ra)>0)
        ra = mean(ra);
    else
        ra = -1;
    end
    if(numel(rb)>0)
        rb = mean(rb);
    else
        rb = -1;
    end
end