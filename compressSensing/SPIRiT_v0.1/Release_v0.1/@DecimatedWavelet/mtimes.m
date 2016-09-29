function res = mtimes(a,b)

if isa(a,'DecimatedWavelet') == 0
    error('In  A.*B only A can be DecimatedWavelet operator');
end

if a.adjoint
    N = size(b,1);
    sizeINI = b{1,1}.sizeINI;
    res = zeros(sizeINI(1), sizeINI(2), N);
else
    N = size(b,3);
    res = cell(N,2);
end

% for every coil
for n=1:N
    if a.adjoint        
        iWT_R = indwt2(b{n,1});
        iWT_I = indwt2(b{n,2});        
        sizeINI = b{n,1}.sizeINI;
        res(:,:,n) = crop(iWT_R,sizeINI(1),sizeINI(2)) + i*crop(iWT_I,sizeINI(1),sizeINI(2));
    else
        % forward transform
        [WT_R, s] = wavedec2(real(b(:,:,n)), a.N, a.wname, a.ExtM); 
        [WT_I, s] = wavedec2(imag(b(:,:,n)), a.N, a.wname, a.ExtM); 
        res{n,1} = WT_R;
        res{n,2} = WT_I;
    end
end
