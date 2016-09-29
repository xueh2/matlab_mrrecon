function res = fft3c_Permute(x)
% x: [COL LIN CHA PAR ...]

s = size(x);

if ( numel(s) == 3 )       
    res = fft3c(x);
end

if ( numel(s) == 4 ) % COL LIN CHA PAR
    x2 = permute(x, [1 2 4 3]);
    res = fft3c(x2);
    res = permute(res, [1 2 4 3]);
end

if ( numel(s) == 5 )
    x2 = permute(x, [1 2 4 3 5]);
    res = fft3c(x2);
    res = permute(res, [1 2 4 3 5]);
end
