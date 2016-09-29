function [A, B, T1, res, Jac, cov, mse] = PerformT1SASHAFitting_nlinfit(TI, y, options)
% perform the t1 sasha fitting using nlinfit
% options.ThreePara: if 1, three para SASHA will be called, if not, two para sasha is called
% options.Robust: 'on' or 'off', whether to perform robust fit

A = max(y(:)) + 1;
T1 = median(TI);

if ( options.ThreePara )
    B = A;
    initialPara = [A B T1];
else
    initialPara = [A T1];
end

if ( options.ThreePara )
    [beta,res,Jac,cov,mse] = nlinfit(TI,y,@modelfuncThreePara,initialPara, options);
else
    [beta,res,Jac,cov,mse] = nlinfit(TI,y,@modelfuncTwoPara,initialPara, options);
end

A = beta(1);

if ( options.ThreePara )
    B = beta(2);
    T1 = beta(3);
else
    B = beta(1);
    T1 = beta(2);
end

end

function yhat = modelfuncTwoPara(b,x)

        if ( abs(b(2)) < eps )
            b(2) = sign(b(2))*eps;
        end

        yhat = b(1) - b(1)*exp(-x./b(2));        
        nans = isnan(yhat(:));
        yhat(find(nans>0)) = 0;
        
        ind = isfinite(yhat);
        yhat(find(ind==0)) = 0;
end

function yhat = modelfuncThreePara(b,x)

        if ( abs(b(3)) < eps )
            b(3) = sign(b(3))*eps;
        end
        
        yhat = b(1) - b(2)*exp(-x./b(3));
        nans = isnan(yhat(:));
        yhat(find(nans>0)) = 0;
        ind = isfinite(yhat);
        yhat(find(ind==0)) = 0;
end
