function [sF, rF, sKi, rKi, sE, rE, sPS, rPS, sVb, rVb] = analyze_res_table(res_table)
% [sF, rF, sKi, rKi, sE, rE, sPS, rPS, sVb, rVb] = analyze_res_table(res_table)

ind = find(res_table.gender=='F');
disp(['Total number :' num2str(numel(res_table.gender))]);
disp(['Female :' num2str(numel(ind))]);
disp(['Age :' num2str(mean(res_table.age)) '+/-' num2str(std(res_table.age))]);
disp('================================================')
N = size(res_table.sf,1);

ind = find(res_table.sf(:)>0);
ind2 = find(res_table.rf(:)>0);
disp(['Stress flow - ' num2str(mean(res_table.sf(ind))) '+/-' num2str(std(res_table.sf(ind))) ' ; rest flow - ' num2str(mean(res_table.rf(ind2))) '+/-' num2str(std(res_table.rf(ind2))) ]);

% plot res table only
ind = find(res_table.sf(:)>0);
ind2 = find(res_table.rf(:)>0);

figure
hold on

plot(res_table.rf(ind2), res_table.rKi_Fermi(ind2), 'b.');
plot(res_table.rf(ind2), res_table.rKi_BTEX(ind2), 'bo');

plot(res_table.sf(ind), res_table.sKi_Fermi(ind), 'r.');
plot(res_table.sf(ind), res_table.sKi_BTEX(ind), 'k+');
plot(res_table.sf(ind), res_table.sKi_MF(ind)*0.9, 'k^');
plot(res_table.sf(ind), res_table.sKi_TwoCompExp(ind), 'ko');

plot([0 7], [0 7], 'k-');
hold off
xlabel('Flow (ml/min/g), BTEX', 'FontSize', 12)
ylabel('Ki (ml/min/g), Fermi', 'FontSize', 12)
box on
grid on
xlim([0 7]);
ylim([0 7]);
legend('Fermi, rest', 'BTEX Ki, rest', 'Fermi, stress', 'BTEX Ki, stress', 'MF Ki, stress', 'TwoCompExp Ki, stress');
title('Flow - Ki plot, BTEX vs. Model-Free', 'FontSize', 12)

for n=1:N
    
    rf = res_table.rf(n, :);
    sf = res_table.sf(n, :);
       
    ind = find(rf>0);
    rF(n) = mean(rf(ind));
    
    ind = find(sf>0);
    sF(n) = mean(sf(ind));
    
    % ----------------------
    rKi_v = res_table.rKi_MF(n, :);
    ind = find(rKi_v>0);
    rKi(n,1) = mean(rKi_v(ind));
    
    rKi_v = res_table.rKi_Fermi(n, :);
    ind = find(rKi_v>0);
    rKi(n,2) = mean(rKi_v(ind));
    
    rKi_v = res_table.rKi_TwoCompExp(n, :);
    ind = find(rKi_v>0);
    rKi(n,3) = mean(rKi_v(ind));
    
    rKi_v = res_table.rKi_BTEX(n, :);
    ind = find(rKi_v>0);
    rKi(n,4) = mean(rKi_v(ind));
    
    % ----------------------
    sKi_v = res_table.sKi_MF(n, :);
    ind = find(sKi_v>0);
    sKi(n,1) = mean(sKi_v(ind));
    
    sKi_v = res_table.sKi_Fermi(n, :);
    ind = find(sKi_v>0);
    sKi(n,2) = mean(sKi_v(ind));
    
    sKi_v = res_table.sKi_TwoCompExp(n, :);
    ind = find(sKi_v>0);
    sKi(n,3) = mean(sKi_v(ind));
    
    sKi_v = res_table.sKi_BTEX(n, :);
    ind = find(sKi_v>0);
    sKi(n,4) = mean(sKi_v(ind));
    
    % ----------------------
    rE_v = res_table.rE(n,:);
    ind = find(rE_v>0 & rE_v<1);
    rE(n) = mean(rE_v(ind));
    
    sE_v = res_table.sE(n,:);
    ind = find(sE_v>0 & sE_v<1);
    sE(n) = mean(sE_v(ind));
    
    % ----------------------
    rPS_v = res_table.rPS(n,:);
    ind = find(rPS_v>0);
    rPS(n) = mean(rPS_v(ind));
    
    sPS_v = res_table.sPS(n,:);
    ind = find(sPS_v>0);
    sPS(n) = mean(sPS_v(ind));
    
    % ----------------------
    rVisf_v = res_table.rVisf(n,:);
    ind = find(rVisf_v>0);
    rVisf(n) = mean(rVisf_v(ind));
    
    sVisf_v = res_table.sVisf(n,:);
    ind = find(sVisf_v>0);
    sVisf(n) = mean(sVisf_v(ind));
    
    % ----------------------
    rVb_v = res_table.rVb(n,:);
    ind = find(rVb_v>0);
    rVb(n) = mean(rVb_v(ind));
    
    sVb_v = res_table.sVb(n,:);
    ind = find(sVb_v>0);
    sVb(n) = mean(sVb_v(ind));
end

disp(['stress flow: ' num2str(mean(sF)) '+/-' num2str(std(sF))]);
disp(['rest flow: ' num2str(mean(rF)) '+/-' num2str(std(rF))]);
disp(['stress E: ' num2str(mean(sE)) '+/-' num2str(std(sE))]);
disp(['rest E: ' num2str(mean(rE)) '+/-' num2str(std(rE))]);
disp(['stress PS: ' num2str(mean(sPS)) '+/-' num2str(std(sPS))]);
disp(['rest PS: ' num2str(mean(rPS)) '+/-' num2str(std(rPS))]);
disp(['stress Visf: ' num2str(mean(sVisf)) '+/-' num2str(std(sVisf))]);
disp(['rest Visf: ' num2str(mean(rVisf)) '+/-' num2str(std(rVisf))]);
disp(['stress Vb: ' num2str(mean(sVb)) '+/-' num2str(std(sVb))]);
disp(['rest Vb: ' num2str(mean(rVb)) '+/-' num2str(std(rVb))]);
disp('================================================')
disp(['stress Ki [MF, Fermi, TwoCompExp, BTEX]: ' num2str(mean(sKi, 1)) '+/-' num2str(std(sKi, [], 1))]);
disp(['rest Ki [MF, Fermi, TwoCompExp, BTEX]: ' num2str(mean(rKi, 1)) '+/-' num2str(std(rKi, [], 1))]);
disp('================================================')
ind = find(res_table.sGd(:)>0);
disp(['stress peak Gd in myo: ' num2str(mean(res_table.sGd(ind))) '+/-' num2str(std(res_table.sGd(ind)))]);
ind = find(res_table.rGd(:)>0);
disp(['rest peak Gd in myo: ' num2str(mean(res_table.rGd(ind))) '+/-' num2str(std(res_table.rGd(ind)))]);
disp('================================================')
ind = find(res_table.sSNR(:)>0);
disp(['stress peak SNR: ' num2str(mean(res_table.sSNR(ind))) '+/-' num2str(std(res_table.sSNR(ind)))]);
ind = find(res_table.rSNR(:)>0);
disp(['rest peak SNR: ' num2str(mean(res_table.rSNR(ind))) '+/-' num2str(std(res_table.rSNR(ind)))]);
disp('================================================')
disp(['stress peak Gd of AIF: ' num2str(mean(res_table.sAifPeakGd)) '+/-' num2str(std(res_table.sAifPeakGd))]);
disp(['rest peak Gd of AIF: ' num2str(mean(res_table.rAifPeakGd)) '+/-' num2str(std(res_table.rAifPeakGd))]);
disp('================================================')
disp(['stress T2* at peak Gd: ' num2str(mean(res_table.sT2S)) '+/-' num2str(std(res_table.sT2S))]);
disp(['rest T2* at peak Gd: ' num2str(mean(res_table.rT2S)) '+/-' num2str(std(res_table.rT2S))]);
disp('================================================')
