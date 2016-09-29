function res = mtimes(a,b)

if isa(a,'UndecimatedWavelet') == 0
    error('In  A.*B only A can be UndecimatedWavelet operator');
end

% for every coil
poolSize = matlabpool('size');

if a.adjoint
    N = size(b,1);
    sizeINI = b{1,1}.sizeINI;
    res = zeros(sizeINI(1), sizeINI(2), N);

    if poolSize > 0 
        parfor n=1:N
            iWT_R = indwt2(b{n,1});
            iWT_I = indwt2(b{n,2});        
            sizeINI = b{n,1}.sizeINI;
            res(:,:,n) = crop(iWT_R,sizeINI(1),sizeINI(2)) + i*crop(iWT_I,sizeINI(1),sizeINI(2));
        end
    else
        for n=1:N
            iWT_R = indwt2(b{n,1});
            iWT_I = indwt2(b{n,2});        
            sizeINI = b{n,1}.sizeINI;
            res(:,:,n) = crop(iWT_R,sizeINI(1),sizeINI(2)) + i*crop(iWT_I,sizeINI(1),sizeINI(2));
        end
    end

else    
    if poolSize > 0
        N = size(b,3);
        res = cell(N,2);
        parfor n=1:N
            WT_R = ndwt2(real(b(:,:,n)), a.N, a.wname, 'mode', a.ExtM);
            res{n,1} = WT_R;
        end
        
        parfor n=1:N
            WT_I = ndwt2(imag(b(:,:,n)), a.N, a.wname, 'mode', a.ExtM);
            res{n,2} = WT_I;
        end        
    else
        N = size(b,3);
        res = cell(N,2);
        for n=1:N
            % forward transform
            WT_R = ndwt2(real(b(:,:,n)), a.N, a.wname, 'mode', a.ExtM);
            WT_I = ndwt2(imag(b(:,:,n)), a.N, a.wname, 'mode', a.ExtM);
            res{n,1} = WT_R;
            res{n,2} = WT_I;
        end
    end    
end
