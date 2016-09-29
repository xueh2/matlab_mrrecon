function res = times(a,b)

if isa(a,'UndecimatedWavelet') == 0
    error('In  A.*B only A can be UndecimatedWavelet operator');
end

res = cell(n,2);

% for every coil
for n=1:size(b,3)
    if a.adjoint        
        iWT_R = indwt2(b{n,1});
        iWT_I = indwt2(b{n,2});        
        sizeINI = b{n,1}.sizeINI;
        res = crop(iWT_R,sizeINI(1),sizeINI(2)) + i*crop(iWT_I,sizeINI(1),sizeINI(2));
    else
        % forward transform
        WT_R = ndwt2(real(b:,:,n), a.N, a.wname, 'mode', a.ExtM);
        WT_I = ndwt2(imag(b:,:,n), a.N, a.wname, 'mode', a.ExtM);
        res{n,1} = WT_R;
        res{n,2} = WT_I;
    end
end


