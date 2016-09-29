function x = cgRegularizedSPIRiTWithLineSearch_GPU(y, GOP, nIter, x0, params, centre, width)
%-----------------------------------------------------------------------
%
% x = cgRegularizedSPIRiTWithLineSearch(y, GOP, nIter, lambda, x0)
%
% implementation of a regularized non linear conjugate gradient reconstruction
%
% The function solves the following problem:
%
% given k-space measurments y, and the SPIRiT operator G
% using the null space concept, the acquired points are y
% the unacquired points are yu
% D is the selection operator for y
% Dc is the selection operator for yu
% W is the wavelet operator -- sparsity transform
% finds the kspace x that minimizes:
%
% Phi(x) = ||(G-I)*Dc'*yu+(G-I)*D'*y||^2 + lambda*|W*F'*(Dc'*yu+D'*y)|_1
%
%
% the optimization method used is non linear conjugate gradient with fast&cheap backtracking
% line-search.
% 
% Hui Xue 2011
%-------------------------------------------------------------------------

tCGSPIRiT = tic;

if nargin<7
    centre = 1024;
    width = 1024;
end

[sx,sy,nCoils] = size(y);

idx_acq = find(abs(y)>0);
idx_acq = gpuArray(idx_acq);

idx_nacq = find(abs(y)==0);
idx_nacq = gpuArray(idx_nacq);

y = gpuArray(y);
x0 = gpuArray(x0);

x = x0;

% data parameters
kernel = getKernel(GOP);
kSize = [size(kernel,1),size(kernel,2)];

[Nfe, Npe, nCoils] = size(y);

N = length(idx_nacq(:));

% store the acquired data
params.GOP = GOP;
params.y = y;
params.idx_acq = idx_acq;
params.idx_nacq = idx_nacq;
params.nsize = [Nfe Npe nCoils];
params.usingGPU = (strcmp(getMethod(GOP),'fft_AllGPU')==1);
params.fftScaling = sqrt(sx*sy);

% line search parameters
maxlsiter = params.lineSearchItnlim ;
gradToll = params.gradToll ;
objToll = params.objToll;
alpha = params.lineSearchAlpha; 
beta = params.lineSearchBeta;
t0 = params.lineSearchT0;
k = 0;

t = 1;

show = params.show;
    
% precompute
tempY = params.y;
calibCostY = params.GOP*tempY;

% initial gradient
g0 = wGradient(x,params);

dx = -g0;

% secant parameters
prevX = x;
sx = x;
sigma0 = 1;
jsecant = 0;
nIterSecantMax = 10;
secantThres = 1e-3;

prevF = inf;
% iterations
while(1)
    
    tic
    
    % backtracking line-search

	% pre-calculate values, such that it would be cheap to compute the objective
	% many times for efficient line-search
%     lsiter = 0;
%     [calibCostYU, calibCostDYU] = preobjective(x, dx, params);
% 	f0 = objective(calibCostY, calibCostYU, calibCostDYU, x, dx, 0, params);
% 	t = t0;
%     [f1, ERRobj, RMSerr] = objective(calibCostY, calibCostYU, calibCostDYU, x, dx, t, params);

%     while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)
% 	% while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:))) & (lsiter<maxlsiter)
% 		lsiter = lsiter + 1;
% 		t = t * beta;
% 		[f1, ERRobj, RMSerr]  =  objective(calibCostY, calibCostYU, calibCostDYU, x, dx, t, params);
%     end

    % fix equation
    % compute A*dx
%     tempY(:) = 0;
%     tempY(idx_nacq) = dx(idx_nacq);
%     v = params.GOP*tempY;
%     t = (g0(:)'*g0(:)) / ( dx(:)'*v(:)+eps ); --> does not work

% 	if lsiter == maxlsiter
% 		disp('Reached max line search,.... not so good... might have a bug in operators. exiting... ');
% 		return;
% 	end
% 
% 	% control the number of line searches by adapting the initial step search
% 	if lsiter > 2
% 		t0 = t0 * beta;
% 	end 
% 	
% 	if lsiter<1
% 		t0 = t0 / beta;
% 	end
% 
% 	x = (x + t*dx);

% --------------------------------------------------

    [calibCostYU, calibCostDYU] = preobjective(x, dx, params);
%  	f0 = objective(calibCostY, calibCostYU, calibCostDYU, x, dx, 0, params);

    % secant line-search
    lsiter = 0;
    sx = x;
    phiPrev = wGradient(x+t0*dx,params);
    phiPrev = phiPrev(:)'*dx(:);
    alpha = -t0;
    deltaD = dx(:)'*dx(:);
    thresValue = alpha'*alpha*deltaD;
    
    if ( params.usingGPU )
        while ( lsiter<nIterSecantMax & gather(thresValue)>secantThres )
            phi = wGradient(sx,params);
            phi = phi(:)'*dx(:);
            % alpha = alpha* abs(phi/(phiPrev-phi));
            alpha = alpha* phi/(phiPrev-phi);
            sx = sx + alpha.*dx;
            phiPrev = phi;
            lsiter = lsiter+1;
            thresValue = alpha'*alpha*deltaD;
        end   
    else
        while ( lsiter<nIterSecantMax & thresValue>secantThres )
            phi = wGradient(sx,params);
            phi = phi(:)'*dx(:);
            % alpha = alpha* abs(phi/(phiPrev-phi));
            alpha = alpha* phi/(phiPrev-phi);
            sx = sx + alpha.*dx;
            phiPrev = phi;
            lsiter = lsiter+1;
            thresValue = alpha'*alpha*deltaD;
        end 
    end
    
	if lsiter == maxlsiter
		disp(sprintf('Time used : %f - Reached max line search,.... not so good... might have a bug in operators. exiting... '), toc);
		return;
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2
		t0 = t0 * beta;
	end 
	
	if lsiter<1
		t0 = t0 / beta;
	end

	% x = (x + t*dx);
    prevX = x;
    x = sx;

    [f1, ERRobj, RMSerr] = objective(calibCostY, calibCostYU, calibCostDYU, x, dx, 0, params);
    % f1 = objective(calibCostY, calibCostYU, calibCostDYU, x, dx, 0, params);
    % ----------------------------------------------------------------
    
    %conjugate gradient calculation
    
	g1 = wGradient(x,params);
	
    % Fletcher - Reeves updates
    bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	
    % Polak - Ribiere updates
    % bk = (g1(:)'*g1(:)-abs(g1(:)'*g0(:)))/(g0(:)'*g0(:)+eps);
    
    % Dixon's formula
    % bk = abs(g1(:)'*g1(:)/(dx(:)'*g0(:)+eps))
    
    % make sure the CG degrades to steepest descent
    % bk = max([0 bk]);
    
    g0 = g1;
	dx =  - g1 + bk.* dx;
	k = k + 1;
	
	%TODO: need to "think" of a "better" stopping criteria ;-)
    if (params.usingGPU)
        f1 = gather(f1);
        dxNorm = gather(sqrt(dx(:)'*dx(:)));
    else
        dxNorm = norm(dx);
    end
    
	if (k > nIter) | (dxNorm<gradToll) | (prevF-f1<objToll)
        disp(sprintf('Time used : %f - Final iteration : %d --- grad: %f, obj reduction: %f', toc, k,dxNorm,prevF-f1));
        if ( prevF-f1<0 )
            x = prevX; % roll back
        end
		break;
    end

    prevF = f1;
    
    if ( mod(k, params.continuationStep) == 0 )
        params.xfmWeight = params.xfmWeight / params.wavWeightRatio;
        params.TVWeight = params.TVWeight / params.wavWeightRatio;
        disp(['xfm weight : ' num2str(params.xfmWeight) ' - ' ' TV weight : ' num2str(params.TVWeight)]);
    end
    
    if (params.usingGPU) 
        disp(sprintf('Time used : %f - Iter %d --- obj: %f, RMS: %f, L-S: %d', toc, k,f1,gather(RMSerr),lsiter));    
    else
        disp(sprintf('Time used : %f - Iter %d --- obj: %f, RMS: %f, L-S: %d', toc, k,f1,RMSerr,lsiter));    
    end
    
    if show
        if ( k == 1 ) 
            h=-1; 
        end
        
        if (params.usingGPU)
            h = plotKSpaceInFigure(gather(x), [1 1 1], centre, width, 1, 1/24, h); 
            title(['cg-NullSpace-GPU * Iter : ' num2str(k) ' xfm: ' num2str(params.xfmWeight) ' - TV: ' num2str(params.TVWeight)]); 
            drawnow;
        else
            h = plotKSpaceInFigure(x, [1 1 1], centre, width, 1, 1/24, h); 
            title(['cg-NullSpace-GPU * Iter : ' num2str(k) ' xfm: ' num2str(params.xfmWeight) ' - TV: ' num2str(params.TVWeight)]); 
            drawnow;
        end
    end        
end

if (params.usingGPU)
    x = gather(x);
end

disp(sprintf('<--------------- Total time used : %f --------------------->', toc(tCGSPIRiT)));

return;


function [calibCostYU, calibCostDYU] = preobjective(yu, dyu, params)

% precalculates transforms to make line search cheap
if ( params.usingGPU )
    tempY = parallel.gpu.GPUArray.zeros(size(params.y));
else
    tempY = zeros(size(params.y));
end

tempY(params.idx_nacq) = yu(params.idx_nacq);
calibCostYU = params.GOP*tempY;

tempY(params.idx_nacq) = dyu(params.idx_nacq);
calibCostDYU = params.GOP*tempY;

% tempY = params.y;
% tempY(params.idx_nacq) = yu(params.idx_nacq);
% Im = ifft2c(tempY);
% sparseCoeff = params.sparseTransform*Im;
% sparseCoeffNormYU = params.sparseTransform.coeffNorm(sparseCoeff);
% sparseCoeffNormYU = sum(sparseCoeffNormYU);

function [res, obj, RMS] = objective(calibCostY, calibCostYU, calibCostDYU, yu, dyu, t, params)
%calculated the objective function

p = params.pNorm;

tempY = params.y;

% tempY(params.idx_nacq) = yu(params.idx_nacq) + t.*dyu(params.idx_nacq);
% obj = params.GOP*tempY;

obj = calibCostYU + t*calibCostDYU + calibCostY;
obj = obj(:)'*obj(:);

if params.xfmWeight  
    tempY(params.idx_nacq) = yu(params.idx_nacq) + t*dyu(params.idx_nacq);
    % Im = ifft2c(tempY);
    Im = params.fftScaling*fftshift(ifft2(ifftshift(tempY)));
    if ( params.usingGPU ) 
        sparseCoeff = gpuArray(params.sparseTransform*gather(Im));
    else
        sparseCoeff = params.sparseTransform*Im;
    end
    [m, sparseCoeffNormYUDYU]= params.sparseTransform.coeffNorm(sparseCoeff);
    XFM = sparseCoeffNormYUDYU;
else
    XFM=0;
end

if params.TVWeight
    tempY(params.idx_nacq) = yu(params.idx_nacq) + t*dyu(params.idx_nacq);
    % Im = ifft2c(tempY);
    Im = params.fftScaling*fftshift(ifft2(ifftshift(tempY)));
    w = params.TV*Im;
    TV = (w(:).*conj(w(:))+params.l1Smooth).^(p/2); 
else
    TV = 0;
end

if params.tikWeight
    tempY(params.idx_nacq) = yu(params.idx_nacq) + t*dyu(params.idx_nacq);
    TIK = tempY(:)'*tempY(:);
else
    TIK = 0;
end

XFM = sum(XFM.*params.xfmWeight(:));
TIK = sum(TIK.*params.tikWeight(:));
TV = sum(TV.*params.TVWeight(:));

res = obj + (XFM) + TIK + TV;
RMS = sqrt(res/sum(abs(params.y(:))>0));

function grad = wGradient(x,params)

gradXFM = 0;
gradTV = 0;
gradTikhonov = 0;

gradObj = gOBJ(x,params);

if params.xfmWeight
    gradXFM = gXFM(x,params);
end

if params.TVWeight
    gradTV = gTV(x,params);
end

if params.tikWeight
    gradTikhonov = gTikhonov(x,params);
end

grad = (gradObj +  params.xfmWeight.*gradXFM + params.tikWeight.*gradTikhonov + params.TVWeight.*gradTV);

function gradObj = gOBJ(yu,params)
% computes the gradient of the data consistency
tempY = params.y;
tempY(params.idx_nacq) = yu(params.idx_nacq);
gradObj = params.GOP'*(params.GOP*tempY);
gradObj = 2*gradObj;

% apply the Dc
gradObj(params.idx_acq) = 0;

function grad = gXFM(yu,params)
% compute gradient of the L1 transform operator
tempY = params.y;
tempY(params.idx_nacq) = yu(params.idx_nacq);

Im = params.fftScaling*fftshift(ifft2(ifftshift(tempY))); 

if ( params.usingGPU )    
    sparseCoeff = gpuArray(params.sparseTransform*gather(Im));
    [sparseCoeffNorm, totalNorm] = params.sparseTransform.coeffNorm(sparseCoeff);
    sparseCoeff = params.sparseTransform.divideByNorm(sparseCoeff, sparseCoeffNorm, params.pNorm, params.l1Smooth);    
    Im = gpuArray(params.sparseTransform'*gather(sparseCoeff));
else
    sparseCoeff = params.sparseTransform*Im;
    [sparseCoeffNorm, totalNorm] = params.sparseTransform.coeffNorm(sparseCoeff);
    sparseCoeff = params.sparseTransform.divideByNorm(sparseCoeff, sparseCoeffNorm, params.pNorm, params.l1Smooth);
    Im = params.sparseTransform'*sparseCoeff;
end

grad = 1/params.fftScaling*fftshift(fft2(ifftshift(Im))); 
grad(params.idx_acq) = 0;

function grad = gTikhonov(yu,params)
% compute gradient of the tikhonov operator
grad = params.y;
grad(params.idx_nacq) = yu(params.idx_nacq);
grad(params.idx_acq) = 0;

function grad = gTV(yu,params)
% compute gradient of TV operator

p = params.pNorm;

tempY = params.y;
tempY(params.idx_nacq) = yu(params.idx_nacq);
% Im = ifft2c(tempY);
Im = params.fftScaling*fftshift(ifft2(ifftshift(tempY)));

Dx = params.TV*Im;

G = p*Dx.*(Dx.*conj(Dx) + params.l1Smooth).^(p/2-1);
grad = params.TV'*G;
% grad = fft2c(grad);
grad = 1/params.fftScaling*fftshift(fft2(ifftshift(grad)));
grad(params.idx_acq) = 0;

