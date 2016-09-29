function [F_L2_BSpline, F_L1, F_L1_BSpline, F_Fermi] = PerformDeconvolution_ParameterSelection(cin, y, r, orderBSpline, numOfInternalControlPoints, sigma, rep)
% perform the monte-carlo to select best parameters
% rep is the number of repitition for every trial
% for every orderBSpline, every number of control points, every sigma (SNR level)
% compute the flow from deconvolution and compare against the best value

%% set up common parameters
N = numel(cin);

if numel(y)~=N
    error('cin and y have different length');
end

F = max(r); % best flow

O = numel(orderBSpline);
M = numel(numOfInternalControlPoints);
S = numel(sigma);

y_ori = y;

y_used = zeros([N rep S]);
for s=1:numel(sigma)
    disp(['--> sigma is ' num2str(sigma(s))]);
    
    for rr=1:rep
        y_used(:, rr, s) = y_ori + sigma(s)*randn(N, 1);
    end
end

F_L2_BSpline = zeros(O, M, S, rep);
F_L1 = zeros(S, rep);
F_L1_BSpline = zeros(O, M, S, rep);
F_Fermi = zeros(O, M, S, rep);

best_lambda_L1 = zeros(S, rep);

maxDistRatio = 0.5;
thres_svd = 1e-4;

for s=1:numel(sigma)
    disp(['--> sigma is ' num2str(sigma(s))]);
    
    for rr=1:rep
        y = y_used(:, rr, s);

        lambda = -1;
        [impulse2, yr2, lambda] = PerformDeconvolution_L1_Fista(cin, y, lambda);
        best_lambda_L1(s, rr) = lambda;
        
        [v, ind] = PerformDeconvolution_findFirstPeak(impulse2, maxDistRatio);       
        F_L1(s, rr) = v;                       
    end
end

for s=1:numel(sigma)
    disp(['--> sigma is ' num2str(sigma(s))]);
    
    % lambda = median( best_lambda_L1(s, :) );
    % lambda = 0.001;
    lambda = 0.01*sigma(s);
    
    for rr=1:rep
        y = y_used(:, rr, s);
               
        for m=1:M
            for o=1:O
                [impulse, yr] = PerformDeconvolution_Tikhonov_BSpline(cin, y, thres_svd, orderBSpline(o), numOfInternalControlPoints(m), lambda);
                [impulse3, yr3] = PerformDeconvolution_L1_Fista_BSpline(cin, y, orderBSpline(o), numOfInternalControlPoints(m), lambda);

                [v, ind] = PerformDeconvolution_findFirstPeak(impulse, maxDistRatio);
                F_L2_BSpline(o, m, s, rr) = v;
                
                [v, ind] = PerformDeconvolution_findFirstPeak(impulse3, maxDistRatio);
                F_L1_BSpline(o, m, s, rr) = v;
                
                [FF, tau, k, td, yr4, impulse4] = PerformDeconvolution_Fermi(cin, y, v);
                F_Fermi(o, m, s, rr) = max(impulse4);
            end
        end        
    end
end

