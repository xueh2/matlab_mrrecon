function filteredData = performRawDataFilter3D(kspace, feFilter, peFilter, parFilter)
% ----------------------------------------------------------
% perform raw data filter
% kspace : kspace data with [Nfe Npe numOfCoils]
% feFilter : filter along the FE direction, [Nfe 1]
% peFilter : filter along the PE direction, [Npe 1]
% parFilter : filter along the PAR direction, [Npar 1]
% ----------------------------------------------------------

filteredData = kspace;        
% numOfPE = size(kspace, 2);
% numOfFE = size(kspace, 1);
% numOfCoils = size(kspace, 3);
% numOfPAR = size(kspace, 4);

COL = size(kspace, 1);
LIN = size(kspace, 2);
CHA = size(kspace, 3);
PAR = numel(parFilter);

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

% for par=1:PAR
%     for pe=1:LIN
%         for fe=1:COL
%             filter3D(fe, pe, par) = feFilter(fe)*peFilter(pe)*parFilter(par);
%         end
%     end
% end

filter3D = reshape(filter3D, [COL LIN 1 PAR]);

for c=1:CHA
    d = kspace(:, :, c, :);                
    filteredData(:, :, c, :) = filter3D.*d;
end
            
% if ( ~isempty(feFilter) )
%     r = 1/sqrt(1/numOfFE * sum(feFilter.*feFilter));
%     for par=1:numOfPAR
%         for pe=1:numOfPE
%             for c=1:numOfCoils
%                 d = kspace(:, pe, c, par);                
%                 filteredData(:, pe, c, par) = r .* feFilter.*d;
%             end
%         end
%     end
% end
% 
% if ( ~isempty(peFilter) )
%     r = 1/sqrt(1/numOfPE * sum(peFilter.*peFilter));
%     for par=1:numOfPAR
%         for fe=1:numOfFE
%             for c=1:numOfCoils
%                 d = filteredData(fe, :, c, par);                
%                 filteredData(fe, :, c, par) = r .* peFilter'.*d;
%             end
%         end
%     end
% end
% 
% if ( ~isempty(parFilter) )
%     r = 1/sqrt(1/numOfPAR * sum(parFilter.*parFilter));
%     for pe=1:numOfPE
%         for fe=1:numOfFE
%             for c=1:numOfCoils
%                 d = filteredData(fe, pe, c, :);                
%                 filteredData(fe, pe, c, :) = r .* parFilter'.*d;
%             end
%         end
%     end
% end
