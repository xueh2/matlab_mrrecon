

function Coef = GRAPPA_SrcDstChannels_Kernel_2D(kspace, acs, option, thresReg)
% function Coef = GRAPPA_SrcDstChannels_Kernel_2D(kspace, acs, option, thresReg)
% kspace : the regularly undersampled kspace [Nfe Npe CHA]
% acs : the acs lines [Nfe Npe CHA], fully sampled acs lines
% option.KernelSize, [k_fe k_pe], kernel size along fe and pe direction 
% option.KernelPattern, kernel lines along PE, e.g. for R4, it could be [-4 0 4 8]
% option.OutPattern, indexes of lines to be generated, e.g. for R4, it could be [1 2 3] 
% option.feRange, the range along FE used for kernel estimation, this field is applied for both kspace and acs
% option.acqPELines, the index of acquired PE lines in kspace, starting from 1
% option.acqAcsLines, the index of acquired acs lines in the coordinate of kspace
%
% Coef : grappa kernel estimated
% the size of Coef is [KernelSize(1)*KernelSize(2)*srcCha length(OutPattern)*dstCha]
% e.g., for a [5 4] kernel, 32 src channels and 4 dst channels, R4, missing lines [1 2 3]
% the coefficient size is [5*4*32 3*4]
% the order of stored coefficients are [fe pe srcCha missingLine dstCha]
% that is, for Coef(1:20,1), the first 5 elements are first column of kernel for srcCha=1
% for Coef(1:20,1), the next 5 elements are second column of kernel for srcCha=1
% for Coef(21:40,1), the first 5 elements are first column of kernel for srcCha=2
% Coef(:, 1:3) are kernels for first dst channel, Coef(:, 4:6) are kernels for the second dst channel etc.

%% apply the fe range

if ( isfield(option, 'feRange') )    
    kspace = kspace(option.feRange(1):option.feRange(2), :,:);
    acs = kspace(option.feRange(1):option.feRange(2), :,:);
end

if ( size(kspace,1) ~= size(acs,1) )    
    feRange = [1 min([size(kspace,1) size(acs,1)])];
    kspace = kspace(feRange(1):feRange(2), :,:);
    acs = kspace(feRange(1):feRange(2), :,:);    
end

%% perpare for the kernel estimation
KernelSize = option.KernelSize;
OutPattern = option.OutPattern;

halfFEKernelSize = floor(KernelSize(1)/2);
feKernelOffset = -halfFEKernelSize:halfFEKernelSize;
if ( numel(feKernelOffset) ~= KernelSize(1) )
    feKernelOffset = -halfFEKernelSize:halfFEKernelSize-1;
end

Start_fe = abs(feKernelOffset(1));
End_fe = feKernelOffset(end);
Pattern_fe = feKernelOffset;

Nfe = size(kspace, 1);
Npe = size(kspace, 2);
srcCha = size(kspace, 3);

% number of acs lines
NpeAcs = size(acs, 2);
% destination channel number
dstCha = size(acs, 3);

% number of grappa coefficients are [KernelSize(1)*KernelSize(2)*srcCha length(OutPattern)*dstCha]
Coef = zeros([KernelSize(1)*KernelSize(2)*srcCha length(OutPattern)*dstCha]);
Coef = complex(Coef, Coef);
CoefSize = prod(KernelSize)*srcCha;

MisSize = length(OutPattern)*dstCha;

% maximal size of data matrix A = [Nfe*NpeAcs KernelSize(1)*KernelSize(2)*srcCha]
% it is possible that some acs lines do not have equations
A = zeros([Nfe*NpeAcs KernelSize(1)*KernelSize(2)*srcCha]);
A = complex(A, A);

% calibration matrix has the maximal size of [Nfe*NpeAcs length(OutPattern)*dstCha]
C = zeros(Nfe*NpeAcs, length(OutPattern)*dstCha);
C = complex(C, C);

% the line offsets for the block patterns to every output lines
lineOffsets = zeros(numel(option.KernelPattern), numel(OutPattern));
for o=1:numel(OutPattern)
    lineOffsets(:,o) = option.KernelPattern - OutPattern(o);
end

% the line offsets for the missing lines
misLineOffsets = OutPattern - OutPattern(1);

%% set up the equation A*Coef = C (acquired lines * kernel = calibration lines)
counter  = 0;
counter2  = 0;

% go through every acs lines and try to put together the equation
for l=1:NpeAcs
    
    % get the block ind of kspace lines
    blockNeeded = lineOffsets(:,1)+option.acqAcsLines(l);
    ind = find(blockNeeded<=0);
    blockNeeded(ind) = blockNeeded(ind) + Npe;
    ind = find(blockNeeded>Npe);
    blockNeeded(ind) = blockNeeded(ind) - Npe;    
    
    % if all needed lines are acquired in kspace, the equation can be listed
    allAcquired = isempty(setdiff(blockNeeded, option.acqPELines))==1;
    
    % if all acs lines are in acs, the equation can be listed
    acsNeeded = l + misLineOffsets;
    ind = find(acsNeeded<=0);
    acsNeeded(ind) = acsNeeded(ind) + Npe;
    ind = find(acsNeeded>Npe);
    acsNeeded(ind) = acsNeeded(ind) - Npe;
    
    allAcsAcquired = isempty(setdiff(acsNeeded, 1:NpeAcs))==1;      
    
    if ( allAcquired & allAcsAcquired )
             
        %% line by line assignment
%         Atemp = kspace(:, blockNeeded(:), :);        
%         for feShift = 1:KernelSize(1)
%             AtempShifted = circshift(Atemp, [-feKernelOffset(feShift) 0 0]); % note the sign of shift
%             for cha=1:srcCha
%                 for pe=1:KernelSize(2)
%                     colInd = feShift + (pe-1)*KernelSize(1) + (cha-1)*KernelSize(1)*KernelSize(2);
%                     A(counter*Nfe+1:(counter+1)*Nfe, colInd) = AtempShifted(:, pe, cha); 
%                 end
%             end
%         end
%              
%         Ctemp = acs(:, acsNeeded(:), :);        
%         for mis=1:numel(OutPattern)
%             for cha=1:dstCha            
%                 colInd = mis + (cha-1)*numel(OutPattern);
%                 C(counter*Nfe+1:(counter+1)*Nfe, colInd) = Ctemp(:, mis, cha);
%             end
%         end
%         
%         counter = counter + 1;
        
        %% point by point assignment
%         for index_fe = Start_fe+1:Nfe-End_fe
%             counter2  = counter2  + 1; 
%             temp = kspace( index_fe+Pattern_fe, blockNeeded, :);
%             A(counter2, :) = temp(:);
%             C(counter2, :) = reshape( acs( index_fe, acsNeeded, :), 1, MisSize);
%         end

        for index_fe = 1:Nfe
            counter2  = counter2  + 1;             
            feInd = index_fe+Pattern_fe;
            feInd = mod(feInd-1, Nfe) + 1;
            temp = kspace( feInd, blockNeeded, :);
            A(counter2, :) = temp(:);
            C(counter2, :) = reshape( acs( index_fe, acsNeeded, :), 1, MisSize);
        end
    end
end

% A = A(1:counter*Nfe, :);
% C = C(1:counter*Nfe, :);
A = A(1:counter2, :);
C = C(1:counter2, :);

Asquare = (A'*A);
AsquareInv = inv_tikhonov_IcePAT(Asquare, thresReg);
Coef = AsquareInv * (A'*C);
% Coef = AsquareInv \ (A'*C);
