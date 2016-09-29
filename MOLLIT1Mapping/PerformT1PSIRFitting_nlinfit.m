function [A, B, T1Star, T1, res, Jac, cov, mse] = PerformT1PSIRFitting_nlinfit(TI, y, options)
% perform the t1 psir fitting using nlinfit

% // A can be initialized as the maximal intensity
% A = *std::max_element(intensities.begin(), intensities.end()) + 1.0;
% 
% // B can be initizlied as the A - min(intensties)
% B = A - *std::min_element(intensities.begin(), intensities.end());
% 
% // From Andreas' MOLLI code, the heuristic rule is that B < m_constrainForB*A
% B = (B<m_constrainForB*A) ? (m_constrainForB*A) : (B+1.0);
% 
% // T1-star is initialized as prepTimes[swapIndex+1]
% // double T1_Star = ( (swapIndex+1) < prepTimes.size() ) ? prepTimes[swapIndex+1] : *prepTimes.end();
% double T1_Star = 0.0;
% 
% if ( swapIndex+1 >= prepTimes.size() )
% {
%     T1_Star = prepTimes[prepTimes.size()-1];
% }
% else
% {
%     if ( FTK_ABS(intensities[swapIndex+1]-intensities[swapIndex]) > FTK_TRUNCATION_DELTA )
%         T1_Star = prepTimes[swapIndex] - (prepTimes[swapIndex+1]-prepTimes[swapIndex])*intensities[swapIndex]/(intensities[swapIndex+1]-intensities[swapIndex]);
%     else
%         T1_Star = 0.5 * (prepTimes[swapIndex] + prepTimes[swapIndex+1]);
% }
% 
% if ( B>FTK_TRUNCATION_DELTA && A>FTK_TRUNCATION_DELTA )
%     T1_Star *= log(B/A);
% else
%     T1_Star = 0;
% 
% // store the initial parameters
% initialPara[0] = A;
% initialPara[1] = B;
% initialPara[2] = T1_Star;

A = max(y(:)) + 1;
B = A - min(y(:));

zeroCrossing = -1;

N = numel(y);
for ind=2:N    
    if ( y(ind-1)*y(ind)<0 )
        zeroCrossing = ind;
        break;
    end    
end

if ( zeroCrossing ~= -1 )    
    value = y(zeroCrossing)-y(zeroCrossing-1);
    if ( abs(value) < eps )
        value = sign(value)*eps;
    end
    x = zeroCrossing-1 - y(zeroCrossing-1)/(value);
    T1Star = interp1(1:N, TI, x);
else
    if ( y(end) <= 0 )
        T1Star = TI(end);
    else
        T1Star = TI(1);
    end
end

initialPara = [A B T1Star];

[beta,res,Jac,cov,mse] = nlinfit(TI,y,@modelfunc,initialPara, options);

A = beta(1);
B = beta(2);
T1Star = beta(3);

if ( abs(A) < eps )
    A = sign(A)*eps;
end
    
T1 = T1Star * (B/A-1);

end

function yhat = modelfunc(b,x)
%     if ( b(3) < eps )
%         yhat = repmat(b(1), [numel(x) 1]);
%     else
        yhat = b(1) - b(2)*exp(-x./b(3));
%     end
end