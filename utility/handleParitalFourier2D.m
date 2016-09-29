function kspaceAfterPF = handleParitalFourier2D(kspace, rangeUsedCOL, rangeUsedLIN)
% handle the particial fourier
% kspace: [COL LIN CHA PHS/REP]
% 'Sym': fill with conjugate symmetry of kspace
% 'Huang': fill with Feng Huang's method
% rangeUsedCOL, LIN are actually acquired kspace range

s = size(kspace);

COL = s(1);
LIN = s(2);
CHA = s(3);
PHS = s(4);

kspaceAfterPF = kspace;

if ( isempty(rangeUsedCOL) || isempty(rangeUsedLIN) )
    return;
end

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

% asymmetric echo
    
for phs=1:PHS
    for col=1:COL           
        reflectedCOL = -(col-1 + COLShift) - COLShift + 1;

        for lin=1:LIN

            if ( col>=rangeUsedCOL(1) & col<=rangeUsedCOL(2) & lin>=rangeUsedLIN(1) & lin<=rangeUsedLIN(2)  )
                continue;
            end

            reflectedLIN = -(lin-1 + LINShift) - LINShift + 1;

            if ( reflectedCOL<=COL & reflectedLIN<=LIN )
                kspaceAfterPF(col, lin, :, phs) = conj(kspace(reflectedCOL, reflectedLIN, :, phs));
            end
        end
    end
end
