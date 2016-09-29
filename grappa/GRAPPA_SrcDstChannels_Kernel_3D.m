

function Coeff = GRAPPA_SrcDstChannels_Kernel_3D(src, dst, option, thresReg)
% function Coef = GRAPPA_SrcDstChannels_Kernel_3D(src, dst, option, thresReg)
% src : the regularly undersampled src [Nfe Npe CHA Npar]
% dst : the dst lines [Nfe Npe CHA Npar], fully sampled dst lines
% option.KernelSize, [k_fe k_pe k_par], kernel size along fe, pe, par direction 
% option.KernelPattern, kernel lines along PE, e.g. for R4, it could be [-4 0 4 8]
% option.KernelPatternPAR, kernel lines along PAR, e.g. for R4, it could be [-4 0 4 8]
% option.OutPattern, indexes of lines to be generated, e.g. for R4, it could be [1 2 3] or [0 1 2 3]
% option.OutPatternPAR, indexes of lines to be generated, e.g. for R4, it could be [1 2 3] or [0 1 2 3]
% option.feRange, the range along FE used for kernel estimation, this field is applied for both src and dst
%
% Coef : grappa kernel estimated
% the size of Coef is [KernelSize(1) KernelSize(2) KernelSize(3) srcCha dstCha length(OutPatternPE) length(OutPatternPAR)]
% e.g., for a [5 4 4] kernel, 32 src channels and 4 dst channels, R4*3, missing lines [1 2 3] for PE and [1 2] for PAR
% the coefficient size is [5 4 4 32 4 3 2]
% the order of stored coefficients are [fe pe par srcCha dstCha ope opar]

%% apply the fe range

if ( isfield(option, 'feRange') )    
    src = src(option.feRange(1):option.feRange(2), :,:,:);
    dst = dst(option.feRange(1):option.feRange(2), :,:,:);
end

if ( size(src,1) ~= size(dst,1) )    
    feRange = [1 min([size(src,1) size(dst,1)])];
    src = src(feRange(1):feRange(2), :,:,:);
    dst = dst(feRange(1):feRange(2), :,:,:);    
end

%% perpare for the kernel estimation
KernelSize = option.KernelSize;

KernelPattern = option.KernelPattern;
KernelPatternPAR = option.KernelPatternPAR;

kPE = numel(KernelPattern);
kPAR = numel(KernelPatternPAR);

OutPatternPE = option.OutPattern;
OutPatternPAR = option.OutPatternPAR;

oPE = numel(OutPatternPE);
oPAR = numel(OutPatternPAR);

halfKfe = floor(KernelSize(1)/2);
KernelSize(1) = 2*halfKfe + 1;

Nfe = size(src, 1);
Npe = size(src, 2);
srcCha = size(src, 3);
Npar = size(src, 4);

% destination channel number
dstCha = size(dst, 3);

% assemble the equations
sFE = halfKfe + 1;
eFE = Nfe - halfKfe;

sPE = abs(KernelPattern(1))+1;
ePE = Npe - KernelPattern(end);

sPAR = abs(KernelPatternPAR(1))+1;
ePAR = Npar - KernelPatternPAR(end);

lenFE = eFE-sFE+1;

% A*K = D
rowA = lenFE * (ePE-sPE+1) * (ePAR-sPAR+1);
colA = srcCha * KernelSize(1) * KernelSize(2) * KernelSize(3);

rowD = rowA;
colD = dstCha * oPE * oPAR;

A = zeros(rowA, colA);
A = complex(A, A);

D = zeros(rowD, colD);
D = complex(D, D);

%% set up the equation A*Coef = D (acquired lines * kernel = calibration lines)
% go through every point

rowInd = 1;
for par=sPAR:ePAR
    for pe=sPE:ePE
        for fe=sFE:eFE
            
            colInd = 1;
            for kpar=1:kPAR
                opar = par + KernelPatternPAR(kpar);

                for kpe=1:kPE
                    ope = pe + KernelPattern(kpe);

%                     for kfe=-halfKfe:halfKfe
% 
%                         ofe = fe + kfe;
% 
%                         A(rowInd, colInd:colInd+srcCha-1) = reshape(src(ofe, ope, :, opar), [1 srcCha]);
%                         colInd = colInd + srcCha;
%                     end

                      A(rowInd, colInd:colInd+KernelSize(1)*srcCha-1) = reshape(src(fe-halfKfe:fe+halfKfe, ope, :, opar), [1 KernelSize(1)*srcCha]);
                      colInd = colInd + KernelSize(1)*srcCha;
                end
            end
           
            colInd = 1;
            for kpar=1:oPAR
                opar = par + OutPatternPAR(kpar);

                for kpe=1:oPE
                    ope = pe + OutPatternPE(kpe);

                    D(rowInd, colInd:colInd+dstCha-1) = reshape(dst(fe, ope, :, opar), [1 dstCha]);
                    colInd = colInd + dstCha;
                end
            end            
            
            rowInd = rowInd+1;
        end
    end
end

Asquare = (A'*A);
AsquareInv = inv_tikhonov_IcePAT(Asquare, thresReg);
K = AsquareInv * (A'*D);

Coeff = zeros(KernelSize(1), KernelSize(2), KernelSize(3), srcCha, dstCha, oPE, oPAR);

for opar=1:oPAR
    for ope=1:oPE
        for dcha=1:dstCha            
            colInd = dcha + (ope-1)*dstCha + (opar-1)*dstCha*oPE;
            
            for kpar=1:kPAR
                for kpe=1:kPE
                    for kfe=1:KernelSize(1)
                        for scha=1:srcCha    
                            % rowInd = scha + (kfe-1)*srcCha + (kpe-1)*KernelSize(1)*srcCha + (kpar-1)*KernelSize(1)*srcCha*kPE;
                            rowInd = kfe + (scha-1)*KernelSize(1) + (kpe-1)*KernelSize(1)*srcCha + (kpar-1)*KernelSize(1)*srcCha*kPE;
                            
                            Coeff(kfe, kpe, kpar, scha, dcha, ope, opar) = K(rowInd, colInd);
                        end
                    end
                end
            end
        end
    end
end


