function filteredData = performRawDataFilter(kspace, feFilter, peFilter)
% ----------------------------------------------------------
% perform raw data filter
% kspace : kspace data with [Nfe Npe numOfCoils]
% feFilter : filter along the FE direction, [Nfe 1]
% peFilter : filter along the PE direction, [Npe 1]
% ----------------------------------------------------------

filteredData = kspace;        
numOfPE = size(kspace, 2);
numOfFE = size(kspace, 1);
numOfCoils = size(kspace, 3);

% remember to scale the filter to make sure they keep the SNR unit

if ( ~isempty(feFilter) )
    r = 1/sqrt(1/numOfFE * sum(feFilter.*feFilter));
    for pe=1:numOfPE
        for c=1:numOfCoils
            d = kspace(:, pe, c);                
            filteredData(:, pe, c) = r .* feFilter.*d;
        end
    end
end

if ( ~isempty(peFilter) )
    r = 1/sqrt(1/numOfPE * sum(peFilter.*peFilter));
    for fe=1:numOfFE
        for c=1:numOfCoils
            d = filteredData(fe, :, c);                
            filteredData(fe, :, c) = r .* peFilter'.*d;
        end
    end
end
