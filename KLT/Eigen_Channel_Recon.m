
function Data = Eigen_Channel_Recon(Data, V, numberOfModesKept)
% perform the eigen channel recon
% Data(:,:,1) corresponds to the smallest eigenvalue
% Data(:,:,end) corresponds to the largest eigenvalue
% numberOfModesKept : number of modes kept for recon; if -1, all mode are used

S0 = size(Data);

if ( length(S0) == 2 )
    numOfFe = S0(1);
    numOfPe = 1;
    numOfCoil = S0(2);

    Data = reshape( Data, numOfFe*numOfPe, numOfCoil );

    if ( numberOfModesKept > 0 )
        V(:, 1:numOfCoil-numberOfModesKept) = 0;
    end

    Data = Data * V';
    Data = reshape( Data, numOfFe, numOfPe, numOfCoil );   
end

if ( length(S0) == 3 )
    numOfFe = S0(1);
    numOfPe = S0(2);
    numOfCoil = S0(3);

    Data = reshape( Data, numOfFe*numOfPe, numOfCoil );

    if ( numberOfModesKept > 0 )
        V(:, 1:numOfCoil-numberOfModesKept) = 0;
    end

    Data = Data * V';
    Data = reshape( Data, numOfFe, numOfPe, numOfCoil );   
end

if ( length(S0) == 4 )
    numOfFe = S0(1);
    numOfPe = S0(2);
    numOfCoil = S0(3);
    numOfFrame = S0(4);
    
    Data = permute(Data, [1 2 4 3]);
    Data = reshape( Data, numOfFe*numOfPe*numOfFrame, numOfCoil );

    if ( numberOfModesKept > 0 )
        V(:, 1:numOfCoil-numberOfModesKept) = 0;
    end

    Data = Data * V';
    Data = reshape( Data, numOfFe, numOfPe, numOfFrame, numOfCoil );   
    Data = permute(Data, [1 2 4 3]);
end

return









