function res = fft3c_ImageDimensions(x)

COL = size(x, 1);
LIN = size(x, 2);
ACQ = size(x, 3);
SLC = size(x, 4);
PAR = size(x, 5);
ECO = size(x, 6);
PHS = size(x, 7);
REP = size(x, 8);
SET = size(x, 9);

res = x;

for set=1:SET
    for rep=1:REP
        for phs=1:PHS
            for eco=1:ECO
                for slc=1:SLC
                    for acq=1:ACQ
                        d = x(:, :, acq, slc, :, eco, phs,  rep, set);
                        res(:, :, acq, slc, :, eco, phs,  rep, set) = fft3c(squeeze(d));
                    end
                end
            end
        end
    end
end