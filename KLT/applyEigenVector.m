
function Data = applyEigenVector(Data, V)
% apply the eigen vector
% Data(:,:,1) corresponds to the smallest eigenvalue
% Data(:,:,end) corresponds to the largest eigenvalue

S0 = size(Data);
numOfFe = S0(1);
numOfPe = S0(2);
numOfCoil = S0(3);

if ( length(S0) == 3 )
    Data = reshape( Data, numOfFe*numOfPe, numOfCoil );
    Data = Data * V;
    Data = reshape( Data, numOfFe, numOfPe, numOfCoil );   
end

if ( length(S0) == 4 )
    Data2 = Data;
    for f=1:S0(4)
        D = reshape( Data(:,:,:,f), numOfFe*numOfPe, numOfCoil );
        D = D * V;
        Data2(:,:,:,f) = reshape( D, numOfFe, numOfPe, numOfCoil );   
    end
    Data = Data2;
end

return









