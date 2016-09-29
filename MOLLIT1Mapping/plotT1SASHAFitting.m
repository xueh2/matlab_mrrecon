function h = plotT1SASHAFitting(TI, y, tx, otherParas, h)
% plot the t1 sasha fitting

beta = [otherParas(1) otherParas(2) tx];

yhat = modelfuncThreePara(beta,TI);

if ( ishandle(h) )
    figure(h);
    hold on
    plot(TI, y, 'b-');
    plot(TI, y, 'b.');
    plot(TI, yhat, 'r-');
    plot(TI, yhat, 'r.');
    hold off
else
    h = figure;
    hold on
    plot(TI, y, 'b-');
    plot(TI, y, 'b.');
    plot(TI, yhat, 'r-');
    plot(TI, yhat, 'r.');
    hold off
end

box on;

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
