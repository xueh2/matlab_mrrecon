function G = spiritKernelSparseMatrixRepresentation(kernel, Nfe, Npe)
% convert the spirit kernels to the sparse matrix representation
% the size of G is Nfe*Npe*NCoil by Nfe*Npe*NCoil
% kernel : the kspace domain spirit kernel
% the 1D vector of kspace is in the storage order of Nfe - Npe - Coil (same as the matlab column-wise order)

ksize = size(kernel);
nSrcCoil = ksize(3);
nDstCoil = ksize(4);
numOfKs = ksize(1)*ksize(2);

gSize = [Nfe*Npe*nDstCoil Nfe*Npe*nSrcCoil];
nzmax = gSize(1)*numOfKs*nSrcCoil;
G = spalloc(gSize(1), gSize(2), nzmax);

% compute the 1D index for the first coil
% the periodic boundary condistion is used
indNfe = zeros(ksize(1), ksize(2), Nfe, Npe);
indNpe = zeros(ksize(1), ksize(2), Nfe, Npe);
ind1DOrder = zeros(ksize(1), ksize(2), Nfe, Npe);

halfKfe = (ksize(1)-1)/2;
halfKpe = (ksize(2)-1)/2;

tic

toc

tic
for pe=1:Npe
    for fe=1:Nfe        
        for kpe=-halfKpe:halfKpe
            for kfe=-halfKfe:halfKfe
                
                feInd = fe + kfe;
                peInd = pe + kpe;
                
                % periodic boundary condition
                if ( feInd < 1 )
                    feInd = Nfe + feInd;
                end
                
                if ( feInd > Nfe )
                    feInd = feInd - Nfe;
                end
            
                if ( peInd < 1 )
                    peInd = Npe + peInd;
                end
                
                if ( peInd > Npe )
                    peInd = peInd - Npe;
                end
                
                indNfe(kfe+halfKfe+1, kpe+halfKpe+1, fe, pe) = feInd;
                indNpe(kfe+halfKfe+1, kpe+halfKpe+1, fe, pe) = peInd;                
                ind1DOrder(kfe+halfKfe+1, kpe+halfKpe+1, fe, pe) = sub2ind([Nfe Npe], feInd, peInd);
            end
        end                
    end
end
toc

% fill the sparse matrix
pointValues = zeros([nzmax 4]);

tic
pointInd = 1;
for dstc=1:nDstCoil    
    dstOffset = (dstc-1)*Nfe*Npe;
    for srcc=1:nSrcCoil

        disp([num2str(dstc) ' - ' num2str(srcc)]);
        
        srcOffset = (srcc-1)*Nfe*Npe;

        aKernel = kernel(:, :, srcc, dstc);
        aKernel = aKernel(:);
        realK = real(aKernel);
        imagK = imag(aKernel);
        
        for pe=1:Npe
            for fe=1:Nfe

                % dst channel
                % rowInd = sub2ind([Nfe Npe nDstCoil], fe, pe, dstc);
                rowInd = fe + (pe-1)*Nfe + dstOffset;
%                 rowInd = indDst;
%                 indDst = indDst + 1;

%                 srcInds = ind1DOrder(:, :, fe, pe);
%                 srcInds = srcInds(:);
                colInds = ind1DOrder(:, :, fe, pe) + srcOffset;
                colInds = colInds(:);
                
                ppInd = pointInd:pointInd+numOfKs-1;
                pointValues(ppInd, 1) = rowInd;
                pointValues(ppInd, 2) = colInds;
                pointValues(ppInd, 3) = realK;
                pointValues(ppInd, 4) = imagK;
                pointInd = pointInd+numOfKs;
                % end
                
                % G(rowInd, colInds(:)) = aKernel;

%                 for kpe=1:ksize(2)                
%                     for kfe=1:ksize(1)                                        
% 
%                         % src channel
%                         colInd = ind1DOrder(kfe, kpe, fe, pe) + offset;                    
%                         G(rowInd, colInd) = kernel(kfe, kpe, srcc, dstc);
%                     end
%                 end
            end
        end
    end
end

G = sparse(pointValues(:,1), pointValues(:,2), pointValues(:,3)+i*pointValues(:,4), gSize(1), gSize(2));

toc
