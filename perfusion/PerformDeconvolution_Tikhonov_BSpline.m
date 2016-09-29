function [r, yr] = PerformDeconvolution_Tikhonov_BSpline(cin, y, thres_svd, orderBSpline, numOfInternalControlPoints, lambda)
% perform the tikhonov deconvolution with BSpline modelling
% thres_svd: the threshold for svd signular value
% numOfInternalControlPoints: number of internal control points, number of
% knots are numOfInternalControlPoints+2*orderBSpline, the clamped BSpline
% form is used
% lambda: tiknohov regularization lamda

%% set up common parameters
N = numel(cin);

if numel(y)~=N
    error('cin and y have different length');
end

% number of control points, make sure the clamped bspline is used
p = numOfInternalControlPoints + orderBSpline;
knots = zeros(1, orderBSpline+p);

stepSize = (N-1) / (numOfInternalControlPoints+1);

knots(orderBSpline+1:orderBSpline+numOfInternalControlPoints) = stepSize*(1:numOfInternalControlPoints);
knots(orderBSpline+numOfInternalControlPoints+1:end) = N-1;

%% set up the matrix A
A = zeros(N, N);
for i=1:N
    for j=i:-1:1
        A(i,j) = cin(i-j+1);
    end
end

%% set up the matrix D
B = bspline_basismatrix(orderBSpline, knots, 0:N-1);
D = A*B;

%% solve the tiknohov regularization
L = get_l (p,1);

[UU,sm,XX] = cgsvd (D,L);
% [reg_corner,rho,eta,reg_param] = l_curve (UU,sm,y, 'Tikh', L, XX, 0, 0);
[reg_corner,rho,eta,reg_param] = l_curve(UU,sm,y, 'Tikh', L, XX);

coeff = tikhonov (UU,sm,XX,y,reg_corner);

r = bspline_deboor(orderBSpline, knots, coeff', 0:N-1);
yr = A*r';
% figure; hold on; plot(cin); plot(y, '.'); plot(yr, 'r+'); hold off

% sigma_ratio = sm(:,1) ./ sm(1,1);
% ind = find(sigma_ratio<thres_svd);
% if ( ~empty(ind) )
%     r2 = tgsvd (UU,sm,XX,b,1:ind(end));
% end

