function grappaCoeff = grappa_2d_coefficients(data, header, acs, headerAcs, kernalIndFE, blockCoord, missingLines, regularization, thresReg, ... 
                            overDetermineRatio, sampledRangeFE, kspaceCenterFE, phaseEncodingUsed, varargin)
%
% grappaCoeff = grappa_2d_coefficients(data, header, acs, headerAcs, kernalIndFE, blockCoord, missingLines, regularization, thresReg, overDetermineRatio, sampledRangeFE, kspaceCenterFE, phaseEncodingUsed, varargin)
%
% This function performs the simple cartesian 2D grappa coefficients computation
% 
% ------------------------ input ------------------------
% data : kspace data, [Nfe Npe numOfCoil]
% header: information header for data
% acs : the auto-callibration lines [Nfe numOfAcsLines numOfCoil]
% headerAcs: information header for acs signal
% kernalIndFE : for 2D grappa kernel, the relative indexes of points used in the fe direction, e.g. for kernel size 5, it can be [-2 -1 0 1 2]
% blockCoord : define the block structure used by the grappa kernel, blockCoord gives the line numbers of all lines in this block
% e.g. for reduction factor 4, block size 4, blockCoord can be [-4 0 4 8]
% missingLines : given the block structure defined in blcokCoord, missingLines shows the lines which grappa coefficients will be computed for
% missingLines(1) must be the first missing line with minimal PE
% regularization : 'svd' or 'tikhonov'
% thresReg : threshold for regulization
% overDetermineRatio : the overdetermine ratio = number of equantions / number of coefficients for grappa 
% kspaceCenterFE : the kspace centre along the FE direction
% paritalFourierFE : 1 if partial fourier is used along the FE direction
% phaseEncodingUsed : the portion of PE lines used, [0 1], 1 means 100% lines used, 0.5 means 50% lines used
% 
% E.g., if we have 12 coil, reduction factor is 4, blcokCoord is [-4 0 4 8], kernalLengthFE is 5, missingLines is [1 2 3]
% then grappaCoeff should have the size 240*36
% the 36 columns store the grappa kernels of every missing lines of every coil
% first column : 1st missing line for coil 1; 2nd : 2nd missing line for coil 1; 3rd : 3rd missing line for coil 1; 4th : 1st missing line for coil 2 etc.
% for every column, the 240 (20 weights per coil * 12coils) weights are stored as 4*12*5 arrays: 
% 1st dimension is line index, corresponding to block lines 1:4
% 2nd dimension is coil, corresponding to 12 coils 1:12
% 3rd dimension is FE shift, corresponds to FE offsets -2:2
% 
% Assumption: all acs lines are acquired missing lines; mod(numOfAcsLines, length(missingLines)) == 0
% ------------------------ output ------------------------
% grappaCoeff : grappa kernel 
% 
% Hui Xue
% Jan 31, 2011
%
% References:   GRAPPA paper
% ========================================================

% retrieve the parameters used in the computation
sizeData = size(data);
Nfe = sizeData(1);
Npe = sizeData(2);
numOfCoil = sizeData(3);
samplingLocationData = header.sampling_location;

samplingLocationAcs = headerAcs.sampling_location;
numOfAcs = size(acs, 2);

minFEUsed = max(header.minFEUsed, sampledRangeFE(1));
maxFEUsed = min(header.maxFEUsed, sampledRangeFE(2));

blockSize = length(blockCoord);
missingLinesSize = length(missingLines);

if ( mod(numOfAcs, missingLinesSize) ~= 0 )
    error('Number of acs lines must include complete missingLine sets ...');
end

kernalLengthFE = length(kernalIndFE);

% for every missingLine set in the asc lines, find the available block pe lines
numOfAcsMissingSet = floor(numOfAcs/missingLinesSize);

%blockInds = zeros(numOfAcsMissingSet, blockSize);
blockInds = [];
blockLineOffsets = missingLines(1) - blockCoord;
acquiredBlocks = [];
for n=1:numOfAcsMissingSet    
    blockStartingIndexes = samplingLocationAcs((n-1)*missingLinesSize+1);
    blockLineInds = blockStartingIndexes - blockLineOffsets; 
    blockLineInds = mod(blockLineInds, Npe);
    blockLineInds(find(blockLineInds==0)) = Npe;
    
    % if all block lines are acquired
    allAcquired = 1;
    for b=1:blockSize
        if ( isempty(find(blockLineInds(b)==samplingLocationData)) )
            allAcquired = 0;
            break;
        end
    end
    
    if ( allAcquired )
        % blockInds(n, :) = blockLineInds;
        blockInds = [blockInds; blockLineInds];
        acquiredBlocks = [acquiredBlocks n];    
    end    
end
totalUsedBlocks = length(acquiredBlocks);

% find the PE line range
minPE = round(Npe*(0.5-phaseEncodingUsed/2));
if ( minPE < 1 )
    minPE = 1;
end

maxPE = round(Npe*(0.5+phaseEncodingUsed/2));
if ( maxPE > Npe )
    maxPE = Npe;
end

blockIndsWithinRange = [];
acquiredBlocksWithinRange = [];
for b=1:totalUsedBlocks    
    blocks =  blockInds(b, :);
    if ( min(blocks) >= minPE & max(blocks)<=maxPE) 
        blockIndsWithinRange = [blockIndsWithinRange; blocks];
        acquiredBlocksWithinRange = [acquiredBlocksWithinRange b];
    end
end

% blockInds = blockIndsWithinRange;
acquiredBlocks = acquiredBlocksWithinRange;
totalUsedBlocks = length(acquiredBlocksWithinRange);

% Two ways of applying Fourier transform is supported.
% 1: 1D IFFT along FE -> GRAPPA Recon -> 1D IFFT along PE
% 2: GRAPPA Recon -> 2D IFFT
% In order to support the first option, 1D IFFT needs to be processed before GRAPPA reconstruction
if header.fft == 1
    for index = 1:numOfCoil
        data(:,:,index) = ifftshift(ifft(ifftshift(data(:,:,index),2), [], 2),2);
    end
end

% range of the FE points used

if ( overDetermineRatio > 0 )
    numOfUnknowns = blockSize*numOfCoil*kernalLengthFE;
    numOfEquations = overDetermineRatio * numOfUnknowns;   
    numOfFEPointsUsed = round(numOfEquations/totalUsedBlocks);
    halfNum = round(numOfFEPointsUsed/2);
    range = kspaceCenterFE-halfNum:kspaceCenterFE+halfNum;
    
    startInd = range(1);
    endInd = range(end);
    if ( startInd < sampledRangeFE(1) )        
        startInd = sampledRangeFE(1);
        endInd = endInd + startInd - range(1);
    end
    
    if ( endInd > sampledRangeFE(2) )
        endInd = sampledRangeFE(2);
        startInd = startInd - (range(end)-sampledRangeFE(2));
    end
    
    minFEUsed = max(startInd, sampledRangeFE(1));
    maxFEUsed = min(endInd, sampledRangeFE(2));
    
    range = minFEUsed:maxFEUsed;
    numOfFEPoints = length(range);    
else
    range = minFEUsed:maxFEUsed;
    numOfFEPoints = length(range);
end

% The core GRAPPA reconstruction code
% A*grappaCoeff = B
% A is the data matrix with size [totalUsedBlocks*numOfFEPoints by blockSize*numOfCoil*kernalLengthFE]
% B is the acs matrix with size [totalUsedBlocks*numOfFEPoints by numOfCoil*missingLinesSize]
% grappaCoeff is grappa coefficient matrix with size [blockSize*numOfCoil*kernalLengthFE by numOfCoil*missingLinesSize]

A = zeros(totalUsedBlocks*numOfFEPoints, blockSize*numOfCoil*kernalLengthFE);
B = zeros(totalUsedBlocks*numOfFEPoints, numOfCoil*missingLinesSize);

% for every blcok
blockCoilA = zeros(Nfe, blockSize*numOfCoil);
blockA = zeros(Nfe, blockSize*numOfCoil*kernalLengthFE);

blockCoilB = zeros(Nfe, numOfCoil*missingLinesSize);

for b=1:totalUsedBlocks
    n = acquiredBlocks(b);
    blockLineInds = blockInds(n, :);
    
    % fill the matrix A
    % for every coil
    for c=1:numOfCoil
        blockCoilA(:, (c-1)*blockSize+1:c*blockSize) = data(:, blockLineInds, c);
    end
    
    % for the 2D grappa kernel, do the circular shift
    for f = 1:kernalLengthFE       
        blockA(:, (f-1)*blockSize*numOfCoil+1:f*blockSize*numOfCoil) = circshift(blockCoilA, [kernalIndFE(f) 0]);
    end    
    A((b-1)*numOfFEPoints+1:b*numOfFEPoints, :) = blockA(range, :);
    
    % fill the matrix B
    for c=1:numOfCoil
        blockCoilB(:, (c-1)*missingLinesSize+1:c*missingLinesSize) = acs(:, (n-1)*missingLinesSize+1:n*missingLinesSize, c);
    end
    B((b-1)*numOfFEPoints+1:b*numOfFEPoints, :) = blockCoilB(range, :);    
end

grappaCoeff = zeros(blockSize*numOfCoil*kernalLengthFE, numOfCoil*missingLinesSize);

if ( hasGPU() )

    A = gpuArray(A);
    Asquare = (A'*A);
    Asquare = gather(Asquare);
    
    AB = gather(A'*gpuArray(B));
    A = gather(A);
else
    Asquare = (A'*A);
    AB = A'*B;
end

if ( isempty(varargin) )
    if ( strcmp(regularization, 'svd') )
        pseudoInv = inv_reg(Asquare, thresReg);
    end

    if ( strcmp(regularization, 'tikhonov') )
        % pseudoInv = inv_tikhonov(Asquare, 0.001);
        pseudoInv = inv_tikhonov_IcePAT(Asquare, thresReg);
    end

    grappaCoeff = pseudoInv * A' * B;
else
    grappaCoeff = choleskyMatrixSolver(regularize_tikhonov_IcePAT(Asquare, thresReg), AB);
end
