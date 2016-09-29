function res = mtimes(a,b)

if isa(a,'RedundantHarrDWT2D') == 0
    error('In  A.*B only A can be RedundantHarrDWT2D operator');
end

if a.adjoint
    % backward transform
    res = Matlab_PerformRedundantHarrDWT2D(double(b), a.N, -1);
else    
    % forward transform
    res = Matlab_PerformRedundantHarrDWT2D(double(b), a.N, 1);
end
