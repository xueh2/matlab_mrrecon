
clear all
close all

%% bolus rest
cd D:\vessel_utilities\mr\MrRecon\perfusion

load bolus_rest

load bolus_stress

N = numel(cin);
A = zeros(N, N);
for i=1:N
    for j=i:-1:1
        A(i,j) = cin(i-j+1);
    end
end
y = A*r;
% save bolus_rest cin y r
% save bolus_stress cin y r

F = max(r)
sigma = [0.01:0.01:0.15];

y_ori = y;
peakSNR = max(y_ori)./sigma;

rep = 20;
orderBSpline = [3 4];
numOfInternalControlPoints = [3 4 5 6 7 8];

O = numel(orderBSpline);
M = numel(numOfInternalControlPoints);
S = numel(sigma);

[F_L2_BSpline, F_L1, F_L1_BSpline, F_Fermi] = PerformDeconvolution_ParameterSelection(cin, y, r, orderBSpline, numOfInternalControlPoints, sigma, rep);

% F_L2_BSpline = squeeze(F_L2_BSpline);
% F_L1 = squeeze(F_L1);
% F_L1_BSpline = squeeze(F_L1_BSpline);

meanF_L2_BSpline = mean(F_L2_BSpline, 4);
meanF_L1 = mean(F_L1, 2);
meanF_L1_BSpline = mean(F_L1_BSpline, 4);
meanF_Fermi = mean(F_Fermi, 4);

diffF_L2_BSpline = F_L2_BSpline - F;
diffF_L1 = F_L1 - F;
diffF_L1_BSpline = F_L1_BSpline - F;
diffF_Fermi = F_Fermi - F;

diffF_L2_BSpline = mean(abs(diffF_L2_BSpline), 4) * 100 / F;
diffF_L1 = mean(abs(diffF_L1), 2) * 100 / F;
diffF_L1_BSpline = mean(abs(diffF_L1_BSpline), 4) * 100 / F;
diffF_Fermi = mean(abs(diffF_Fermi), 4) * 100 / F;

% find the best bspline knots number

start = find(peakSNR>8)

[bestOrder_L2_ind, bestOrder_L2, bestNumKnots_L2_ind, bestNumKnots_L2] = PerfSimuData_FindBestOrder_NumberKnots(diffF_L2_BSpline, orderBSpline, numOfInternalControlPoints, start)
[bestOrder_L1_ind, bestOrder_L1, bestNumKnots_L1_ind, bestNumKnots_L1] = PerfSimuData_FindBestOrder_NumberKnots(diffF_L1_BSpline, orderBSpline, numOfInternalControlPoints, start)
[bestOrder_Fermi_ind, bestOrder_Fermi, bestNumKnots_Fermi_ind, bestNumKnots_Fermi] = PerfSimuData_FindBestOrder_NumberKnots(diffF_Fermi, orderBSpline, numOfInternalControlPoints, start)

% for the best knots number, plot the error
figure;
hold on
plot(peakSNR, squeeze(diffF_L2_BSpline(bestOrder_L2_ind, bestNumKnots_L2_ind, :)), 'r-');
plot(peakSNR, squeeze(diffF_L1_BSpline(bestOrder_L1_ind, bestNumKnots_L1_ind, :)), 'b--');
plot(peakSNR, squeeze(diffF_Fermi(bestOrder_Fermi_ind, bestNumKnots_Fermi_ind, :)), 'k--');
hold off
xlabel('SNR');
ylabel('error in percentage');
legend('L2 BSpline', 'L1 BSpline', 'Fermi');

%% plots
lambda = 0.001;

y_ori = y;

cin = cin(28:end);
y = y(28:end);

rep = 256;
peakSNR = [1.5:1:30];
sigma = max(y_ori)./peakSNR;
S = numel(sigma);

diff = zeros(S, 4, rep);
for s=1:S
    s
    
    lambda = 0.005 * sigma(s);
    useAbsoluteLambda = true;
    
    for rr=1:rep
        y = y_ori + sigma(s)*randn(numel(y), 1);

        [impulse, yr] = PerformDeconvolution_Tikhonov_BSpline(cin, y, 0.0001, orderBSpline(bestOrder_L2_ind), numOfInternalControlPoints(bestNumKnots_L2_ind), lambda);
        [impulse2, yr2] = PerformDeconvolution_L1_Fista(cin, y, lambda, useAbsoluteLambda);
        [impulse3, yr3] = PerformDeconvolution_L1_Fista_BSpline(cin, y, orderBSpline(bestOrder_L1_ind), numOfInternalControlPoints(bestNumKnots_L1_ind), lambda, useAbsoluteLambda);

        [v, ind] = PerformDeconvolution_findFirstPeak(impulse3, 1);
        [FF, tau, k, td, yr4, impulse4] = PerformDeconvolution_Fermi(cin, y, v);

        [v, ind] = PerformDeconvolution_findFirstPeak(impulse, 0.9);
        F1 = v;

        [v, ind] = PerformDeconvolution_findFirstPeak(impulse2, 0.9);
        F2 = v;

        [v, ind] = PerformDeconvolution_findFirstPeak(impulse3, 0.9);
        F3 = v;

        F4 = max(impulse4);
        
        diff(s, 1, rr) = F1-F;
        diff(s, 2, rr) = F2-F;
        diff(s, 3, rr) = F3-F;
        diff(s, 4, rr) = F4-F;
        
    end
end

diffM = mean(diff, 3);
figure; plot(peakSNR, diffM);legend('L2', 'L1', 'L1 BSpline', 'Fermi');
figure; plot(peakSNR, abs(diffM));legend('L2', 'L1', 'L1 BSpline', 'Fermi');
figure; plot(peakSNR, abs(diffM)*100/F);legend('L2', 'L1', 'L1 BSpline', 'Fermi');

figure; hold on; plot(cin, '.--'); plot(y, 'r+'); plot(y_ori, 'k+'); hold off

figure; hold on; 
plot(cin, '.--'); 
plot(y, 'r+'); 
plot(yr, '-x'); 
plot(yr2, 'g-x'); 
plot(yr3, 'k-x'); 
plot(yr4, 'M-x'); 
hold off
legend('cin', 'y', 'bspline', 'L1 norm', 'L1 norm, bspline', 'Fermi');

figure; 
hold on; 
plot(impulse, 'b-'); 
plot(impulse2, 'g-x'); 
plot(impulse3, 'k-x'); 
plot(impulse4, 'M-x'); 
hold off

[v, ind] = PerformDeconvolution_findFirstPeak(impulse, 0.9);
F1 = v;

[v, ind] = PerformDeconvolution_findFirstPeak(impulse2, 0.9);
F2 = v;

[v, ind] = PerformDeconvolution_findFirstPeak(impulse3, 0.9);
F3 = v;

F4 = max(impulse4);

F1
F2
F3
F4
