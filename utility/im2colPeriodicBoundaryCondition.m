function D = im2colPeriodicBoundaryCondition(A, kSize, maxKSize)
% A is M*N data matrix [Nfe Npe]
% kSize is [kfe kpe], kernel size
% D is the data matrix [M*N kfe*kpe]
% for every data point, the kfe*kpe block is picked and reformat as a vector along col
% the boundary condition is the periodic

Nfe = size(A,1);
Npe = size(A,2);

kfe = kSize(1);
kpe = kSize(2);

hkfe = floor(kfe/2);
hkpe = floor(kfe/2);

maxhkfe = floor(maxKSize(1)/2);
maxhkpe = floor(maxKSize(2)/2);

D = zeros(Nfe*Npe, kfe*kpe);

ptInd = 1;
% for pe=1:Npe
for pe=maxhkfe+1:Npe-maxhkfe
    
    peRange = [pe-hkpe:pe+hkpe];

%     ind = find(peRange<1);
%     if ( ~isempty(ind) )
%         peRange(ind) = Npe - abs(peRange(ind));
%     end
% 
%     ind = find(peRange>Npe);
%     if ( ~isempty(ind) )
%         peRange(ind) = peRange(ind) - Npe;
%     end
            
    % for fe=1:Nfe
    for fe=maxhkpe+1:Nfe-maxhkpe
        
        feRange = [fe-hkfe:fe+hkfe];
       
%         ind = find(feRange<1);
%         if ( ~isempty(ind) )
%             feRange(ind) = Nfe - abs(feRange(ind));
%         end
%         
%         ind = find(feRange>Nfe);
%         if ( ~isempty(ind) )
%             feRange(ind) = feRange(ind) - Nfe;
%         end
        
        [X, Y] = meshgrid(peRange, feRange);
        
        % indA = sub2ind(size(A), Y(:), X(:));
        indA = (X(:)-1)*Nfe + Y(:);
        
        indUsed = find(indA>=1 & indA<=Nfe*Npe);
        indA2 = indA(indUsed);
        
        b = A(indA2(:));
        D(ptInd, indUsed(:)) = b(:);
        
        % D(ptInd, :) = b(:);
        ptInd = ptInd + 1;
    end
end
