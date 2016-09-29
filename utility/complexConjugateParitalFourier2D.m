function kspaceConj = complexConjugateParitalFourier2D(kspace, rangeUsedCOL, rangeUsedLIN)
% handle the particial fourier
% kspace: [Col Line Cha Acq Slice Partition Echo Phase Rep Set Seg]
% rangeUsedCOL, LIN are actually acquired kspace range

COL = size(kspace, 1);
LIN = size(kspace, 2);
CHA = size(kspace, 3);
ACQ = size(kspace, 4);
SLC = size(kspace, 5);
PAR = size(kspace, 6);
ECO = size(kspace, 7);
PHS = size(kspace, 8);
REP = size(kspace, 9);
SET = size(kspace, 10);

kspaceConj = zeros(size(kspace));

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
for set=1:SET
    for rep=1:REP
        for phs=1:PHS
            for eco=1:ECO
                for par=1:PAR
                    for slc=1:SLC
                        for acq=1:ACQ
                            %for cha=1:CHA
                                for col=1:COL           
                                    reflectedCOL = -(col-1 + COLShift) - COLShift + 1;

                                    for lin=1:LIN

                                        reflectedLIN = -(lin-1 + LINShift) - LINShift + 1;

                                        if ( reflectedCOL<=COL & reflectedLIN<=LIN )
                                            kspaceConj(col, lin,:,:,:,:,:,:,:,:) = conj(kspace(reflectedCOL, reflectedLIN,:,:,:,:,:,:,:,:));
                                        end
                                    end
                                end
                            %end
                        end
                    end
                end
            end
        end
    end
end
