
% blood_F3p5_PSg0_noInputforEandP
% blood_F1_PSg0_noInputforEandP

blood_F1_PSg1_noInputforEandP_FermiPaper
blood_F3p5_PSg1_noInputforEandP_FermiPaper

RestCin_blood_F0p6_PSg1_noInputforEandP_FermiPaper
RestCin_blood_F1_PSg1_noInputforEandP_FermiPaper
RestCin_blood_F3p5_PSg1_noInputforEandP_FermiPaper

RestCin_blood_F1_PSg1_Vp0p04_noInputforEandP_FermiPaper
RestCin_blood_F3p5_PSg1_Vp0p04_noInputforEandP_FermiPaper

RestCin_blood_F1_PSg4_Vp0p01_noInputforEandP_FermiPaper
RestCin_blood_F3p5_PSg4_Vp0p01_noInputforEandP_FermiPaper

RestCin_blood_F1_PSg4_Vp0p09_noInputforEandP_FermiPaper
RestCin_blood_F4_PSg4_Vp0p09_noInputforEandP_FermiPaper

RestCin_blood_F1_PSg2_Vp0p09_noInputforEandP_FermiPaper

RestCin_blood_F1_PSg4_Vp0p04_Visf0p6_noInputforEandP_FermiPaper
RestCin_blood_F2_PSg4_Vp0p04_Visf0p6_noInputforEandP_FermiPaper
RestCin_blood_F3_PSg4_Vp0p04_Visf0p6_noInputforEandP_FermiPaper
RestCin_blood_F4_PSg4_Vp0p04_Visf0p6_noInputforEandP_FermiPaper

deltaT = Cout_v__1(2)-Cout_v__1(1)

figure;
hold on
plot(Cin_vascular__1, Cin_vascular);
plot(Cout_v__1, Cout_v, 'k');
plot(Cout_p__1, Cout_p, 'r');
plot(Cout_e__1, Cout_e, 'g');
hold off

sum(Cin_vascular)
sum(Cout_v)
sum(Cout_p)
sum(Cout_e)

cin = Cin_vascular';
y = Cout_v';

cin = cin(1:35);
y = y(1:35);

thres_svd = 1e-4;
orderBSpline = 4;
numOfInternalControlPoints = 4;
numOfInternalControlPoints_L1 = 4;
lambda = 0.0001;

% y2 = y;
% 
% y = y2 + 0.01*randn(numel(y), 1);

[impulse, yr] = PerformDeconvolution_Tikhonov_BSpline(cin, y, thres_svd, orderBSpline, numOfInternalControlPoints, lambda);
[impulse2, yr2] = PerformDeconvolution_L1_Fista(cin, y, lambda);
[impulse3, yr3] = PerformDeconvolution_L1_Fista_BSpline(cin, y, orderBSpline, numOfInternalControlPoints_L1, lambda);

[v, ind] = PerformDeconvolution_findFirstPeak(impulse3, 1);
[FF, tau, k, td, yr4, impulse4] = PerformDeconvolution_Fermi(cin, y, v);
            
figure; hold on; plot(cin, '.--'); plot(y, 'r+'); hold off

figure; hold on; 
plot(cin, '.--'); 
plot(y, 'r+'); 
plot(yr, '-x'); 
plot(yr2, 'g-x'); 
plot(yr3, 'k-x');
plot(yr4, 'm-x');
hold off
legend('cin', 'y', 'bspline', 'L1 norm', 'L1 norm, bspline', 'Fermi');

figure; 
hold on; 
plot(impulse, 'b-'); 
plot(impulse2, 'g-x'); 
plot(impulse3, 'k-x'); 
plot(impulse4, 'm-x'); 
hold off

[v, ind] = PerformDeconvolution_findFirstPeak(impulse, 0.9);
F = max(impulse)

[v, ind] = PerformDeconvolution_findFirstPeak(impulse2, 0.9);
F2 = max(impulse2)

[v, ind] = PerformDeconvolution_findFirstPeak(impulse3, 0.9);
F3 = max(impulse3)

F4 = max(impulse4)

F
F2
F3
F4
         
ratio = 60/deltaT
F = F*ratio
F2 = F2*ratio
F3 = F3*ratio
F4 = F4*ratio

%% bolus rest
cd D:\vessel_utilities\mr\MrRecon\perfusion

load bolus_rest

load bolus_stress

% save bolus_rest cin y r
% save bolus_stress cin y r

F = max(r)

sigma = [0.01:0.005:0.15];
rep = 20;
orderBSpline = 4;
numOfInternalControlPoints = [3 4 5 6 7 8];

[F_L2_BSpline, F_L1, F_L1_BSpline] = PerformDeconvolution_ParameterSelection(cin, y, r, orderBSpline, numOfInternalControlPoints, sigma, rep);

diffF_L2_BSpline = F_L2_BSpline - F;
diffF_L1 = F_L1 - F;
diffF_L1_BSpline = F_L1_BSpline - F;

diffF_L2_BSpline = mean(abs(diffF_L2_BSpline), 4) * 100 / F;
diffF_L1 = mean(abs(diffF_L1), 2) * 100 / F;
diffF_L1_BSpline = mean(abs(diffF_L1_BSpline), 4) * 100 / F;

peakSNR = max(y_ori)./sigma;
k = 6;
figure;
hold on
plot(peakSNR, squeeze(diffF_L2_BSpline(1, k, :))-F, 'r');
plot(peakSNR, squeeze(diffF_L1(:, :)), 'g');
plot(peakSNR, squeeze(diffF_L1_BSpline(1, k, :)), 'k');
hold off

error_L2 = abs(L1_B-F)*100/F;
error_L1_B = abs(L2-F)*100/F;

for m=1:M
    plot(sigma, abs(L1_B(m,:)-F)*100/F, 'r-');
    plot(sigma, abs(L2(m,:)-F)*100/F, 'b--');
end

figure;
hold on
for m=1:M
    plot(sigma, abs(L1_B(m,:)-F)*100/F, 'r-');
    plot(sigma, abs(L2(m,:)-F)*100/F, 'b--');
end
hold off
title(num2str(numOfInternalControlPoints(m)));

y2 = y;

rep = 20;
diff = zeros(3, numel(sigma), rep);
for p=1:numel(sigma)
    p
    
    for rr=1:rep
        y = y2 + sigma(p)*randn(numel(y), 1);

        [impulse, yr] = PerformDeconvolution_Tikhonov_BSpline(cin, y, thres_svd, orderBSpline, numOfInternalControlPoints, lambda);
        % [impulse2, yr2] = PerformDeconvolution_L1_Fista(cin, y, lambda);
        [impulse3, yr3] = PerformDeconvolution_L1_Fista_BSpline(cin, y, orderBSpline, numOfInternalControlPoints_L1, lambda);

        diff(1, p, rep) = max(impulse) - max(r);
        % diff(2, p, rep) = max(impulse2) - max(r);
        diff(3, p, rep) = max(impulse3) - max(r);
    end
end

diff2 = abs(diff);
diff2 = 100* squeeze(mean(diff2, 3)) ./ max(r);

figure;
hold on
plot(sigma, diff2(1, :), 'r');
plot(sigma, diff2(2, :), 'g');
plot(sigma, diff2(3, :), 'k');
hold off

figure; hold on; plot(cin, '.--'); plot(y_ori, 'r+'); plot(y, 'r+'); hold off

figure; hold on; 
plot(cin, '.--'); 
plot(y, 'r+'); 
plot(y2, 'y+'); 
plot(yr, '-x'); 
plot(yr2, 'g-x'); 
plot(yr3, 'k-x'); 
hold off
legend('cin', 'y', 'bspline', 'L1 norm', 'L1 norm, bspline');

figure; 
hold on; 
plot(impulse, 'b-'); 
plot(impulse2, 'g-x'); 
plot(impulse3, 'k-x'); 
hold off

cin2 = cin(1:42);
y2 = y(1:42);
[F, tau, k, td, yr, r] = PerformDeconvolution_Fermi(cin2, y2, 0.015, 25.0, 1, 0);
figure; hold on; plot(cin2, '.--'); plot(y2, 'r+'); plot(yr, 'rx'); hold off

figure; plot(r);

%% save data for MMID4
