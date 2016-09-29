function filteredData = performKSpaceFilter2D(kspace, feFilter, peFilter)
% ----------------------------------------------------------
% perform kspace filter
% kspace : kspace data with [COL, LIN, ACQ, SLC, PAR, ECO, PHS, REP, SET]
% feFilter : filter along the FE direction, [Nfe 1]
% peFilter : filter along the PE direction, [Npe 1]
% ----------------------------------------------------------

filteredData = kspace;        

if ( isempty(feFilter) & isempty(peFilter) )
    return;
end

COL = size(kspace, 1);
LIN = size(kspace, 2);
ACQ = size(kspace, 3);
SLC = size(kspace, 4);
PAR = size(kspace, 5);
ECO = size(kspace, 6);
PHS = size(kspace, 7);
REP = size(kspace, 8);
SET = size(kspace, 9);

% remember to scale the filter to make sure they keep the SNR unit

if ( isempty(feFilter) )
    feFilter = ones(COL, 1);
end
r = 1/sqrt(1/COL * sum(feFilter.*feFilter));
feFilter = feFilter * r;

if ( isempty(peFilter) )
    peFilter = ones(LIN, 1);
end
r = 1/sqrt(1/LIN * sum(peFilter.*peFilter));
peFilter = peFilter * r;

filter2D = feFilter * peFilter';

for set=1:SET
    for rep=1:REP
        for phs=1:PHS
            for eco=1:ECO
                for par=1:PAR
                    for slc=1:SLC
                        for acq=1:ACQ
                            d = kspace(:, :, acq, slc, par, eco, phs,  rep, set);
                            filteredData(:, :, acq, slc, par, eco, phs,  rep, set) = filter2D.*d;
                        end
                    end
                end
            end
        end
    end
end    
