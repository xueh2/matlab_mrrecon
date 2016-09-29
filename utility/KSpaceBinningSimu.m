function pmuTime = KSpaceBinningSimu(acqDuration, RR, TRofOnePhs, E1, accelFactor, alternating)
% simulate the binning acquisition
% pmuTime : [E1/accelFactor PHS], the pmuTime for every readout line

PHS = floor(acqDuration/TRofOnePhs);
LIN = floor(E1/accelFactor);

TR = TRofOnePhs/LIN; % TR for every readout line

pmuTime = zeros(E1, PHS);

for phs=1:PHS
    
    startE1 = mod(phs-1, accelFactor) + 1;
    lin = startE1:accelFactor:E1;
    
    if ( alternating )
        if ( mod(phs, 2) == 0 )
            ind = 0;
            for e1=1:numel(lin)       
                acqTime = (ind+1)*TR + (phs-1)*TRofOnePhs;

                n = floor(acqTime/RR);

                pmuTime(E1-lin(ind+1)+1, phs) = acqTime - n*RR;

                ind = ind + 1;
            end
        else
            ind = 0;
            for e1=1:numel(lin)      
                acqTime = (ind+1)*TR + (phs-1)*TRofOnePhs;

                n = floor(acqTime/RR);

                pmuTime(lin(ind+1), phs) = acqTime - n*RR;

                ind = ind + 1;
            end
        end
    else
        ind = 0;
        for e1=1:numel(lin)      
            acqTime = (ind+1)*TR + (phs-1)*TRofOnePhs;

            n = floor(acqTime/RR);

            pmuTime(lin(ind+1), phs) = acqTime - n*RR;

            ind = ind + 1;
        end
    end
end
