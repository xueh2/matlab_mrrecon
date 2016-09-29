function res = mtimes(a,b)

if isa(a,'RedundantHarrDWT3D') == 0
    error('In  A.*B only A can be RedundantHarrDWT3D operator');
end

% b is a [COL LIN CHA REP] array
s = size(b);
if a.adjoint
    % backward transform
    if ( numel(s) > 3 )
%         res = zeros(s(1)/2, s(2)/2, s(3), s(4)/2);
%         for c=1:size(b, 3)
%             resCurr = Matlab_PerformRedundantHarrDWT3D(double(squeeze(b(:,:,c,:))), a.N, -1);
%             res(:,:,c,:) = resCurr;
%         end
        res = Matlab_PerformRedundantHarrDWT3D(single(b), a.N, -1);
    else        
        res = Matlab_PerformRedundantHarrDWT3D(single(b), a.N, -1);
    end
else    
    % forward transform
    if ( numel(s) > 3 )
        %res = zeros(2*s(1), 2*s(2), s(3), 2*s(4));
%         for c=1:size(b, 3)
%             resCurr = Matlab_PerformRedundantHarrDWT3D(double(squeeze(b(:,:,c,:))), a.N, 1);
%             res(:,:,c,:) = resCurr;
%         end
        res = Matlab_PerformRedundantHarrDWT3D(single(b), a.N, 1);
    else
        b = reshape(b, [s(1) s(2) 1 s(3)]);
        res = Matlab_PerformRedundantHarrDWT3D(single(b), a.N, 1);
    end
end
