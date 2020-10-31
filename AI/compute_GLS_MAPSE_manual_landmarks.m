function res = compute_GLS_MAPSE_manual_landmarks(A, B, C)
% res = compute_GLS_MAPSE_manual_landmarks(A, B, C)

    rA = compute_for_one(A);
    rB = compute_for_one(B);
    rC = compute_for_one(C);

    res = [rA; rB; rC];

end

function r = compute_for_one(A)

    r = [];
    pids = [];
    
    N = size(A, 1);
    
    for k=1:N
        pids = [pids; A{k,2}];
    end
    
    p = unique(pids);

    for k=1:numel(p)
        pid = p(k);
        
        ind = find(pids==pid);
        
        phases = zeros(numel(ind), 1);
        series_str = [];
        for s=1:numel(ind)
            phases(s) = A{ind(s), 4};
            %series_str{s} = A{ind(s), 3};
            
            has_str = 0;
            for tt=1:numel(series_str)
                if(strcmp(A{ind(s), 3}, series_str{tt}))
                    has_str = 1;
                    break;
                end
            end
            if(~has_str)
                series_str = [series_str; {A{ind(s), 3}}];
            end
        end
        
        if(numel(ind)<2)
            continue
        end
        
        num_seris = numel(series_str);
        if(num_seris==1)        
            pts = zeros(3, 2, numel(ind));
            for s=1:numel(ind)
                pts(1, :, s) = A{ind(s), 5};
                pts(2, :, s) = A{ind(s), 6};
                pts(3, :, s) = A{ind(s), 7};
            end

            [GLS, MAPSE, LV_Lengh_PT1, LV_Lengh_PT2, L] = compute_GLS(pts, 0);        

            r = [r; {A{ind(1),1}}, A{ind(1),2}, {A{ind(1),3}}, phases, max(GLS), max(MAPSE)];
        else
            GLS_s = zeros(numel(series_str), 1);
            MAPSE_s = zeros(numel(series_str), 1);
            
            for tt=1:numel(series_str)
                phases = [];
                ind_series = [];
                for s=1:numel(ind)
                    if(strcmp(A{ind(s), 3}, series_str{tt}))
                        phases(s) = A{ind(s), 4};
                        ind_series = [ind_series; ind(s)];
                    end
                end
                
                pts = zeros(3, 2, numel(ind_series));
                for s=1:numel(ind_series)
                    pts(1, :, s) = A{ind_series(s), 5};
                    pts(2, :, s) = A{ind_series(s), 6};
                    pts(3, :, s) = A{ind_series(s), 7};
                end

                [GLS, MAPSE, LV_Lengh_PT1, LV_Lengh_PT2, L] = compute_GLS(pts, 0);
                
                GLS_s(tt) = max(GLS);
                MAPSE_s(tt) = max(MAPSE);
            end
            
            r = [r; {A{ind(1),1}}, A{ind(1),2}, {A{ind(1),3}}, phases, mean(GLS_s), mean(MAPSE_s)];
        end
    end
end
