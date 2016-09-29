function kspaceAfterPF = handleParitalFourier(kspace, rangeUsedCOL, rangeUsedLIN, rangeUsedPAR)
% handle the particial fourier
% kspace: [COL LIN CHA PAR]
% 'Sym': fill with conjugate symmetry of kspace
% 'Huang': fill with Feng Huang's method
% rangeUsedCOL, LIN, PAR are actually acquired kspace range

s = size(kspace);

COL = s(1);
LIN = s(2);
CHA = s(3);

if ( numel(s) >= 4 )
    PAR = s(4);
else
    PAR = 1;
end

kspaceAfterPF = kspace;

if ( mod(COL, 2) == 1 )
    COLShift = -(COL-1)/2;
else
    COLShift = -COL/2;
end

if ( mod(LIN, 2) == 1 )
    LINShift = -(LIN-1)/2;
else
    LINShift = -LIN/2;
end

if ( mod(PAR, 2) == 1 )
    PARShift = -(PAR-1)/2;
else
    PARShift = -PAR/2;
end

if ( isempty(rangeUsedCOL) )
    rangeUsedCOL = [1 COL];
end

if ( isempty(rangeUsedLIN) )
    rangeUsedLIN = [1 LIN];
end

if ( isempty(rangeUsedPAR) )
    rangeUsedPAR = [1 PAR];
end

% asymmetric echo

if ( PAR > 1 )
    
    for par=1:PAR
        par
        reflectedPAR = -(par-1 + PARShift) - PARShift + 1;
        
        for col=1:COL           
            reflectedCOL = -(col-1 + COLShift) - COLShift + 1;

            for lin=1:LIN

                if ( col>=rangeUsedCOL(1) & col<=rangeUsedCOL(2) & lin>=rangeUsedLIN(1) & lin<=rangeUsedLIN(2) & par>=rangeUsedPAR(1) & par<=rangeUsedPAR(2)  )
                    continue;
                end

                reflectedLIN = -(lin-1 + LINShift) - LINShift + 1;

                if ( reflectedCOL<=COL & reflectedLIN<=LIN & reflectedPAR<=PAR )
                    kspaceAfterPF(col, lin, :, par) = conj(kspace(reflectedCOL, reflectedLIN, :, reflectedPAR));
                end
            end
        end
    end
else
          
    for col=1:COL           
        reflectedCOL = -(col-1 + COLShift) - COLShift + 1;

        for lin=1:LIN

            if ( col>=rangeUsedCOL(1) & col<=rangeUsedCOL(2) & lin>=rangeUsedLIN(1) & lin<=rangeUsedLIN(2) )
                continue;
            end

            reflectedLIN = -(lin-1 + LINShift) - LINShift + 1;

            if ( reflectedCOL<=COL & reflectedLIN<=LIN )
                kspaceAfterPF(col, lin, :) = conj(kspace(reflectedCOL, reflectedLIN, :));
            end
        end
    end
    
end

