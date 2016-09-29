function res = mtimes(a,b)

if isa(a,'MotionCorrection2DPlusT') == 0
    error('In A.*B only A can be MotionCorrection2DPlusT operator');
end

% b is a [COL LIN CHA REP] array
s = size(b);
if a.adjoint
    % backward transform
    if ( numel(s) > 3 )
        res = zeros(s(1), s(2), s(3), s(4));
        for c=1:s(3)
            resCurr = PerformComplexImageWarpping(double(squeeze(b(:,:,c,:))), a.header, a.dxInv, a.dyInv, a.keyFrame, a.interpolator);
            res(:,:,c,:) = resCurr;
        end
    else
        res = PerformComplexImageWarpping(double(squeeze(b)), a.header, a.dxInv, a.dyInv, a.keyFrame, a.interpolator);
    end
else    
    % forward transform
    if ( numel(s) > 3 )
        res = zeros(s(1), s(2), s(3), s(4));
        for c=1:s(3)
            resCurr = PerformComplexImageWarpping(double(squeeze(b(:,:,c,:))), a.header, a.dx, a.dy, a.keyFrame, a.interpolator);
            res(:,:,c,:) = resCurr;
        end
    else
        res = PerformComplexImageWarpping(double(squeeze(b)), a.header, a.dx, a.dy, a.keyFrame, a.interpolator);
    end
end
