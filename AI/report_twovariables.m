function report_twovariables(v1, v2, desp, is_intra)

ind = find(v1>0 & v2>0);

r2 = compute_rsquared(v1(ind), v2(ind));

[h, p_test] = ttest2(v1(ind), v2(ind));

figure; 
[means,diffs,meanDiff,CR,linFit] = BlandAltman2(v1(ind)', v2(ind)', 2, '*');
ylim([-5.2 5.2])
xlim([0.0 30.0])
title(desp)
xlabel('mean')
ylabel('diff')

disp(desp)
disp(['mean diff : ' num2str(mean(v1(ind) - v2(ind))) ' +/- ' num2str(std(v1(ind) - v2(ind)))])
disp(['r2 = ' num2str(r2)]);
disp(['p = ' num2str(p_test)]);

a = v1(ind);
b = v2(ind);

SD = sqrt(sum((a-b).^2)/(2*numel(a)));
2*SD;
M = mean([a; b]);
CV = 100 * SD/M;
disp(['CV = ' num2str(CV)]);

if(is_intra)
    % one way ANOVA
    [p,tbl,stats] = anova1([a, b]);
    SD2 = sqrt(tbl{4, 2}/(numel(a)-1));

    % sample size
    nout = sampsizepwr('t',[M SD], M+1, 0.90);

    M = [a b];
    [r, LB, UB, F, df1, df2, p] = ICC(M, '1-k', 0.90, 0.05);

    SEM = SD2 * sqrt(1-r);
    MDC_90 = 1.65 * sqrt(2) * SEM;
   
    disp(['number of samples = ' num2str(nout)]);
    disp(['MDC_90 = ' num2str(MDC_90)]);
end
