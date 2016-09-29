function plotSumSquareReconPDF(n, M, A)
% plot the sum-square-recon probability density function
% x - M/sigma, given the A/sigma
% M is the sum-square value
% A is the real intensity value in the absence of noise

% assume sigma = 1

NM = numel(M);
NA = numel(A);

p = zeros(NM, NA);

for a=1:NA
    if(A(a)==0)
        if(n==1)
            pa = M .* exp(-M.^2/2);
        else
            pa = 0;
        end
    else
        pa = A(a)*((M/A(a)).^n).*exp(-( A(a)*A(a) + M.*M ) ./ 2).* besseli(n-1, A(a)*M);
    end
    p(:,a) = pa(:);
end

figure; hold on
plot(M, p);
box on
hold off
xlabel('SoS SNR M/sigma')
ylabel('Probability')
title(['PDF, Sum Square Recon A = ' num2str(A)])

p2 = zeros(NA, 1);

for tt=1:NA
    Q = M(:) .* p(:,tt);
    p2(tt) = sum(Q(:))*(M(2)-M(1));
end

figure; hold on
plot(A, p2');
box on
hold off
xlabel('True SNR')
ylabel('SoS estimated SNR')
title(['SNR estimation using SoS'])

    function pa = pm(m)
        if(a==0)
            if(n==1)
                pa = m .* exp(-m.^2/2);
            else
                pa = 0;
            end
        else
            pa = a*((m/a).^n).*exp(-( a*a + m*m ) ./ 2).* besseli(n-1, a*m);
        end
        pa = pa .* m;
    end

end
