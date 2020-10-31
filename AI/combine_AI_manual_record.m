function record_AI = combine_AI_manual_record(record_AI, record_P1, record_P2)
% combine records

    N = size(record_AI, 1);

    ptA = record_AI.pt_A;
    ptB = record_AI.pt_B;

    N1 = size(record_P1, 1);
    N2 = size(record_P2, 1);

    p1 = zeros(N1, 1);
    p2 = zeros(N2, 1);

    for k=1:N1
        p1(k) = record_P1{k, 2};
    end

    for k=1:N2
        p2(k) = record_P2{k, 2};
    end

    GLS_A_P1 = zeros(N,1);
    GLS_B_P1 = zeros(N,1);

    GLS_A_P2 = zeros(N,1);
    GLS_B_P2 = zeros(N,1);

    MAPSE_A_P1 = zeros(N,1);
    MAPSE_B_P1 = zeros(N,1);

    MAPSE_A_P2 = zeros(N,1);
    MAPSE_B_P2 = zeros(N,1);

    for k=1:N
        pid_A = record_AI.pt_A(k);
        pid_B = record_AI.pt_B(k);

        [GLS_A_P1(k), MAPSE_A_P1(k)] = compute_G_M(p1, pid_A, record_P1, 'P1');
        [GLS_A_P2(k), MAPSE_A_P2(k)] = compute_G_M(p2, pid_A, record_P2, 'P2');
        
        [GLS_B_P1(k), MAPSE_B_P1(k)] = compute_G_M(p1, pid_B, record_P1, 'P1');
        [GLS_B_P2(k), MAPSE_B_P2(k)] = compute_G_M(p2, pid_B, record_P2, 'P2');

    end
    
    record_AI.GLS_A_P1 = GLS_A_P1;
    record_AI.GLS_A_P2 = GLS_A_P2;
    record_AI.GLS_B_P1 = GLS_B_P1;
    record_AI.GLS_B_P2 = GLS_B_P2;
    record_AI.MAPSE_A_P1 = MAPSE_A_P1;
    record_AI.MAPSE_A_P2 = MAPSE_A_P2;
    record_AI.MAPSE_B_P1 = MAPSE_B_P1;
    record_AI.MAPSE_B_P2 = MAPSE_B_P2;
end

function [G, M] = compute_G_M(p1, pid_A, record_P1, oper_str)
    ind1 = find(p1==pid_A);
    if(isempty(ind1))
        G=-1;
        M=-1;
        return;
    end
    G = [];
    M = [];
    for t=1:numel(ind1)
        G(t) = record_P1{ind1(t), 5};
        M(t) = record_P1{ind1(t), 6};
    end
    ind = find(M>0);
    if((max(M(ind))-min(M(ind)))>0.5 * mean(M(ind)))
        disp([oper_str]);
        record_P1(ind1, :)
    end
    
    M = mean(M(ind));
    
    ind = find(G>0);
    G = mean(G(ind));    
    if(M>30)
        G=-1;
        M=-1;
    end
end