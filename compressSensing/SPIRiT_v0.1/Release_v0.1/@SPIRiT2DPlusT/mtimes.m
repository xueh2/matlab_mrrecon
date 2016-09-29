function res = mtimes(a,x)
% This method performs the SPIRiT operator on an input 2D+t kspace [COL LIN CHA REP]
rep = size(x, 4);

if ( strcmp(a.method, 'fft_AllGPU') )
    res = parallel.gpu.GPUArray.zeros(size(x), 'single');
    
    if ( a.numOfRep == rep )
        if a.adjoint
            for r=1:rep        
                res(:,:,:,r) = a.SpiritGOPs{r}'*x(:,:,:,r); 
            end
        else
            for r=1:rep        
                res(:,:,:,r) = a.SpiritGOPs{r}*x(:,:,:,r); 
            end
        end
    else
        if a.adjoint
            for r=1:rep        
                res(:,:,:,r) = a.SpiritGOPs{1}'*x(:,:,:,r); 
            end
        else
            for r=1:rep        
                res(:,:,:,r) = a.SpiritGOPs{1}*x(:,:,:,r); 
            end
        end        
    end
else
    res = zeros(size(x));
    if ( a.numOfRep == rep )
        if a.adjoint
            for r=1:rep        
                res(:,:,:,r) = a.SpiritGOPs{r}'*x(:,:,:,r); 
            end
        else
            for r=1:rep        
                res(:,:,:,r) = a.SpiritGOPs{r}*x(:,:,:,r); 
            end
        end
    else
        if a.adjoint
            for r=1:rep        
                res(:,:,:,r) = a.SpiritGOPs{1}'*x(:,:,:,r); 
            end
        else
            for r=1:rep        
                res(:,:,:,r) = a.SpiritGOPs{1}*x(:,:,:,r); 
            end
        end        
    end
end
