

function kspaceDst = GRAPPA_SrcDstChannels_Recon_2D(kspace, Coef, option)
% function recon = GRAPPA_SrcDstChannels_Recon_2D(kspace, Coef, option)
% kspace : the regularly undersampled kspace [Nfe Npe CHA]
% Coef : the estiamted grappa kernel using Coef = GRAPPA_SrcDstChannels_Kernel_2D(kspace, acs, option, thresReg)
% option.KernelSize, [k_fe k_pe], kernel size along fe and pe direction 
% option.KernelPattern, kernel lines along PE, e.g. for R4, it could be [-4 0 4 8]
% option.OutPattern, indexes of lines to be generated, e.g. for R4, it could be [1 2 3] 
% option.acqPELines, the index of acquired PE lines in kspace, starting from 1
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
% ---------------------------------------------------------------------------------------------------------
% Note: this program only works for number of destimation channel is equal to number of input channels
% ---------------------------------------------------------------------------------------------------------

s_0 = size(kspace);

KernelSize = option.KernelSize;
Pattern_pe = option.KernelPattern;
OutPattern = option.OutPattern;

halfFEKernelSize = floor(KernelSize(1)/2);
feKernelOffset = -halfFEKernelSize:halfFEKernelSize;
if ( numel(feKernelOffset) ~= KernelSize(1) )
    feKernelOffset = -halfFEKernelSize:halfFEKernelSize-1;
end

Start_fe = abs(feKernelOffset(1));
End_fe = feKernelOffset(end);
Pattern_fe = feKernelOffset;

OutSize = length(OutPattern);
SamplingPattern = double(option.acqPELines) ;
dstCha = floor(size(Coef, 2)/OutSize);

Data_Acq = zeros( s_0(1)*length(SamplingPattern), prod(KernelSize)*s_0(3), 'single' );
Data_Acq = complex( Data_Acq, Data_Acq );

counter  = 0;
for index_pe = 1:length(SamplingPattern)
    k_pe = mod(SamplingPattern(index_pe) + Pattern_pe - 1, s_0(2) ) + 1;
    for index_fe = 1:s_0(1)
        counter  = counter  + 1 ;
        k_fe = mod(index_fe+Pattern_fe-1, s_0(1) ) + 1;
        temp = kspace( k_fe, k_pe, :) ;
        Data_Acq(counter, :) = temp(:);
    end
end

% Reconstruction
Data_Mis = Data_Acq * Coef;

% Reshape to fill k-space
kspaceDst = zeros(s_0(1), s_0(2), dstCha);
if ( dstCha == s_0(3) )
    kspaceDst = kspace;
end

counter = 0;
for index_pe = 1:length(SamplingPattern)
    k_pe = mod( SamplingPattern(index_pe) + OutPattern -1, s_0(2)) + 1;
    for index_fe = 1:s_0(1)    
        counter  = counter  + 1 ;        
        kspaceDst( index_fe, k_pe, :) = reshape( Data_Mis(counter, :), OutSize, dstCha );
    end
end

