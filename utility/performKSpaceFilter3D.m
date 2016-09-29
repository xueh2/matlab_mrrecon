function filteredData = performKSpaceFilter3D(kspace, feFilter, peFilter, parFilter)
% ----------------------------------------------------------
% perform kspace filter
% kspace : kspace data with [COL, LIN, ACQ, SLC, PAR, ECO, PHS, REP, SET]
% feFilter : filter along the FE direction, [Nfe 1]
% peFilter : filter along the PE direction, [Npe 1]
% ----------------------------------------------------------

filteredData = kspace;        

if ( isempty(feFilter) & isempty(peFilter) & isempty(parFilter) )
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

if ( isempty(parFilter) )
    parFilter = ones(PAR, 1);
end
r = 1/sqrt(1/PAR * sum(parFilter.*parFilter));
parFilter = parFilter * r;

filter2D = feFilter * peFilter';
filter3D = repmat(filter2D, [1 1 PAR]);

for par=1:PAR
    filter3D(:,:,par) = filter3D(:,:,par) * parFilter(par);
end

for set=1:SET
    for rep=1:REP
        for phs=1:PHS
            for eco=1:ECO
                for slc=1:SLC
                    for acq=1:ACQ
                        d = kspace(:, :, acq, slc, :, eco, phs,  rep, set);
                        filteredData(:, :, acq, slc, :, eco, phs,  rep, set) = filter3D.* squeeze(d);
                    end
                end
            end
        end
    end
end    
