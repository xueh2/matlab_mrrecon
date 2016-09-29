function fullKspace = grappa_2d_recon(data, header, acs, headerAcs, grappaCoeff, kernalIndFE, blockCoord, missingLines, firstDataLineInBlcok, useACSLines, varargin)
%
% fullKspace = grappa_2d_recon(data, header, acs, headerAcs, grappaCoeff, kernalIndFE, blockCoord, missingLines, varargin)
%
% This function performs the simple cartesian 2D grappa reconstruction
% 
% ------------------------ input ------------------------
% data : kspace data, [Nfe Npe numOfCoil]
% header: information header for data
% acs : the auto-callibration lines [Nfe numOfAcsLines numOfCoil]
% headerAcs: information header for acs signal
% grappaCoeff : grappa coefficents computed using grappa_2d_coefficients
% kernalIndFE : for 2D grappa kernel, the relative indexes of points used in the fe direction, e.g. for kernel size 5, it can be [-2 -1 0 1 2]
% blockCoord : define the block structure used by the grappa kernel, blockCoord gives the line numbers of all lines in this block
% e.g. for reduction factor 4, block size 4, blockCoord can be [-4 0 4 8]
% missingLines : given the block structure defined in blcokCoord, missingLines shows the lines which grappa coefficients will be computed for
% missingLines(1) must be the first missing line with minimal PE
% firstDataLineInBlcok: the index of first data line in data block, e.g. for a block [-4 0 4 8], firstDataLineInBlcok=2 means first data line is matched to 
% line 0 in the block
% useACSLines : if 1, the asc lines are filled into the fullKspace and replace the estimated lines
% 
% E.g., if we have 12 coil, reduction factor is 4, blcokCoord is [-4 0 4 8], kernalIndFE is [-2 -1 0 1 2], missingLines is [1 2 3]
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
% fullKspace : full kspace 
% 
% Hui Xue
% Feb 07, 2011
%
% References:   GRAPPA paper
% ========================================================

% The GRAPPA reconstruction
% A*grappaCoeff = B
% A is the data matrix with size [totalUsedBlocks*numOfFEPoints by blockSize*numOfCoil*kernalLengthFE]
% B is the acs matrix with size [totalUsedBlocks*numOfFEPoints by numOfCoil*missingLinesSize]
% grappaCoeff is grappa coefficient matrix with size [blockSize*numOfCoil*kernalLengthFE by numOfCoil*missingLinesSize]

% retrieve the parameters used in the computation
sizeData = size(data);
Nfe = sizeData(1);
Npe = sizeData(2);
numOfCoil = sizeData(3);
samplingLocationData = header.sampling_location;
lenOfAcquiredData = length(samplingLocationData);

samplingLocationAcs = headerAcs.sampling_location;
numOfAcs = size(acs, 2);

blockSize = length(blockCoord);
missingLinesSize = length(missingLines);

if ( mod(numOfAcs, missingLinesSize) ~= 0 )
    error('Number of acs lines must include complete missingLine sets ...');
end

kernalLengthFE = length(kernalIndFE);

blockLineIndsInData = zeros(lenOfAcquiredData, blockSize);
blockOffset = blockCoord - blockCoord(firstDataLineInBlcok);

missingLineIndsInData = zeros(lenOfAcquiredData, missingLinesSize);
missingLineOffset = missingLines - blockCoord(firstDataLineInBlcok);
for b=1:lenOfAcquiredData
    dataLineInd = samplingLocationData(b);
    blockLineInd = dataLineInd + blockOffset;    
    blockLineIndsInData(b, :) = blockLineInd;    
    
    missingLineIndsInData(b, :) = dataLineInd + missingLineOffset;
end
blockLineIndsInData = mod(blockLineIndsInData, Npe);
blockLineIndsInData(find(blockLineIndsInData==0)) = Npe;

missingLineIndsInData = mod(missingLineIndsInData, Npe);
missingLineIndsInData(find(missingLineIndsInData==0)) = Npe;

% Two ways of applying Fourier transform is supported.
% 1: 1D IFFT along FE -> GRAPPA Recon -> 1D IFFT along PE
% 2: GRAPPA Recon -> 2D IFFT
% In order to support the first option, 1D IFFT needs to be processed before GRAPPA reconstruction
if header.fft == 1
    for index = 1:numOfCoil
        data(:,:,index) = ifftshift(ifft(ifftshift(data(:,:,index),2), [], 2),2);
    end
end

% data matrix
% A*grappaCoeff = B
% B are estimated missing lines for one block for all coils
% this matrix multiplication needs to be repeated lenOfAcquiredData times to estimated all missing lines

A = zeros(Nfe, blockSize*numOfCoil*kernalLengthFE);
B = zeros(Nfe, numOfCoil*missingLinesSize);

% initial the fullKspace
fullKspace = data;
blockCoilA = zeros(Nfe, blockSize*numOfCoil);
for b=1:lenOfAcquiredData
    
    blockLineInd = blockLineIndsInData(b, :);
    missingLineInd = missingLineIndsInData(b, :);
    
    % fill matrix A
    for c=1:numOfCoil
        blockCoilA(:, (c-1)*blockSize+1:c*blockSize) = data(:, blockLineInd, c);
    end
    
    % for the 2D grappa kernel, do the circular shift
    for f = 1:kernalLengthFE       
        A(:, (f-1)*blockSize*numOfCoil+1:f*blockSize*numOfCoil) = circshift(blockCoilA, kernalIndFE(f));
    end    
    
    % compute matrix B
     B = A * grappaCoeff;
     
    % put the estimated lines into fullKspace 
    for c=1:numOfCoil
        fullKspace(:, missingLineInd, c) = B(:, (c-1)*missingLinesSize+1:c*missingLinesSize);
    end
end

% fill the data lines into fullKspace
fullKspace(:, samplingLocationData, :) = data(:, samplingLocationData, :);

% if required, fill the asc lines into the fullKspace
if ( useACSLines )
    for a=1:numOfAcs
        fullKspace(:, samplingLocationAcs(a), :) = acs(:, a, :);
    end
end
