function x = cgRegularizedSPIRiTWithLineSearch_2DPlusT_WithoutNullSpace_GPU(y, GOP, nIter, x0, params, centre, width)
%-----------------------------------------------------------------------
%
% x = cgRegularizedSPIRiTWithLineSearch_2DPlusT_WithoutNullSpace_GPU(y, GOP, nIter, lambda, x0)
%
% implementation of a regularized non linear conjugate gradient reconstruction
%
% The function solves the following problem:
%
% given k-space measurments y, and the SPIRiT operator G
% do not use the null space concept, the acquired points are y
% the unacquired points are yu
% D is the selection operator for y
% Dc is the selection operator for yu
% W is the wavelet operator -- sparsity transform
% finds the kspace x that minimizes:
%
% Phi(x) = ||(G-I)*x||^2 + lambda*|W*F'*x|_1 + lamda2 *||D*x-y||_2 
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

[Nfe, Npe, nCoils, nReps] = size(y);

idx_acq = find(abs(y)>0);
idx_acq = gpuArray(idx_acq);

idx_nacq = find(abs(y)==0);
idx_nacq = gpuArray(idx_nacq);

y = gpuArray(y);
x0 = gpuArray(x0);

x = x0; clear x0

% store the acquired data
params.GOP = GOP; clear GOP
params.y = y;
params.idx_acq = idx_acq;
params.idx_nacq = idx_nacq;
params.nsize = [Nfe Npe nCoils nReps];
params.usingGPU = 1;
params.fftScaling = sqrt(Nfe*Npe);

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
% tempY = params.y;
% calibCostY = params.GOP*tempY;

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
secantRatio = 2;

prevF = inf;
% iterations
while(1)
    
    tic
    
    % [calibCostYU, calibCostDYU] = preobjective(x, dx, params);

    % secant line-search
    lsiter = 0;
    sx = x;
    phiPrev = wGradient(x+t0*dx,params);
    
    norm(phiPrev(:))
    
    phiPrev = phiPrev(:)'*dx(:);
    alpha = -t0;
    deltaD = dx(:)'*dx(:);
    
    deltaD
    
    thresValue = alpha'*alpha*deltaD;
    
    thresValue
    
    prevThresValue = abs(gather(thresValue));
    
    if ( params.usingGPU )
%         while ( lsiter<nIterSecantMax & gather(thresValue)>secantThres )
%             phi = wGradient(sx,params);
%             phi = phi(:)'*dx(:);
%             % alpha = alpha* abs(phi/(phiPrev-phi));
%             alpha = alpha* phi/(phiPrev-phi);
%             sx = sx + alpha.*dx;
%             phiPrev = phi;
%             lsiter = lsiter+1;
%             thresValue = alpha'*alpha*deltaD;
%         end

        while ( lsiter<nIterSecantMax & abs(gather(thresValue))>secantThres & (abs(gather(thresValue))<=secantRatio*prevThresValue) )
            
            phi = wGradient(sx,params);
            
            norm(phi(:))
            
            phi = phi(:)'*dx(:);
            % alpha = alpha* abs(phi/(phiPrev-phi));
            alpha = alpha* phi/(phiPrev-phi);
            
            phiPrev = phi;
            lsiter = lsiter+1;
            prevThresValue = abs(gather(thresValue));
            thresValue = alpha'*alpha*deltaD;
            
            thresValue
            
            if ( abs(gather(thresValue)) <= secantRatio*prevThresValue  )
                sx = sx + alpha.*dx;
                
                norm(sx(:))
                
            end
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

    [f1, ERRobj, RMSerr] = objective(x, dx, 0, params);
    % [f1, ERRobj, RMSerr] = objective(calibCostY, calibCostYU, calibCostDYU, x, dx, 0, params);
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
        disp(['xfm weight : ' num2str(params.xfmWeight) ' - ' ' TV weight : ' num2str(params.TVWeight) ' - ' ' data weight : ' num2str(params.dataWeight)]);
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
            h = plotKSpaceInFigure(gather(x(:,:,:,1)), [1 1 1], centre, width, 1, 1/24, h); 
            title(['cg-NoNullSpace-GPU * Iter : ' num2str(k) ' xfm: ' num2str(params.xfmWeight) ' - TV: ' num2str(params.TVWeight) ' - Data: ' num2str(params.dataWeight)]); 
            drawnow;
        else
            h = plotKSpaceInFigure(x(:,:,:,1), [1 1 1], centre, width, 1, 1/24, h); 
            title(['cg-NoNullSpace-GPU * Iter : ' num2str(k) ' xfm: ' num2str(params.xfmWeight) ' - TV: ' num2str(params.TVWeight) ' - Data: ' num2str(params.dataWeight)]); 
            drawnow;
        end        
    end        
end

if (params.usingGPU)
    x = gather(x);
end

disp(sprintf('<--------------- Total time used : %f --------------------->', toc(tCGSPIRiT)));

return;

function [res, obj, RMS] = objective(yu, dyu, t, params)
%calculated the objective function

p = params.pNorm;

tempY = yu + t.*dyu;

obj = params.GOP*tempY;
obj = obj(:)'*obj(:);

if params.xfmWeight  
    Im = params.fftScaling*fftshift(ifft2(ifftshift(tempY)));
    if ( params.usingGPU ) 
        % sparseCoeff = gpuArray(params.sparseTransform*gather(Im));
        sparseCoeff = params.sparseTransform*gather(Im);
    else
        sparseCoeff = params.sparseTransform*Im;
    end
    
    if ( params.temporalScalingFactor ~= 1 )
        % sparseCoeff(:,:,:, params.nsize(4)+1:end) = params.temporalScalingFactor*sparseCoeff(:,:,:, params.nsize(4)+1:end);
        sparseCoeff = params.sparseTransform.scaleTemporalCoeff(sparseCoeff, params.temporalScalingFactor);
    end
    
    [m, sparseCoeffNormYUDYU]= params.sparseTransform.coeffNorm(sparseCoeff);
    XFM = sparseCoeffNormYUDYU;
else
    XFM=0;
end

if params.TVWeight
    Im = params.fftScaling*fftshift(ifft2(ifftshift(tempY)));
    w = params.TV*Im;
    TV = (w(:).*conj(w(:))+params.l1Smooth).^(p/2); 
else
    TV = 0;
end

if params.tikWeight
    TIK = tempY(:)'*tempY(:);
else
    TIK = 0;
end

if params.dataWeight
    DFED = tempY - params.y;
    DFED(params.idx_nacq) = 0;
    DFED = DFED(:)'*DFED(:);
else
    DFED=0;
end

XFM = sum(XFM.*params.xfmWeight(:));
TIK = sum(TIK.*params.tikWeight(:));
TV = sum(TV.*params.TVWeight(:));
DFED = sum(DFED.*params.dataWeight(:));

res = obj + (XFM) + TIK + TV + DFED;
RMS = sqrt(res/sum(abs(params.y(:))>0));

function grad = wGradient(x,params)

gradXFM = 0;
gradTV = 0;
gradTikhonov = 0;
gradDFED = 0;

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

if params.dataWeight
    gradDFED = gDFED(x,params);
end

grad = (gradObj +  params.xfmWeight.*gradXFM + params.tikWeight.*gradTikhonov ... 
    + params.TVWeight.*gradTV + params.dataWeight.*gradDFED);

function gradObj = gOBJ(yu,params)
% computes the gradient of the data consistency
gradObj = params.GOP'*(params.GOP*yu);
gradObj = 2*double(gradObj);

function gradDFED = gDFED(yu,params);
% data fidelity
% tempY = zeros(size(params.y));
% tempY(params.idx_acq) = yu(params.idx_acq) - params.y(params.idx_acq);
tempY = yu - params.y;
tempY(params.idx_nacq) = 0;
gradDFED = 2*tempY;

function grad = gXFM(yu,params)
% compute gradient of the L1 transform operator
tempY = yu;

Im = params.fftScaling*fftshift(ifft2(ifftshift(tempY))); 

if ( params.usingGPU )    
    % sparseCoeff = gpuArray(params.sparseTransform*gather(Im));
    sparseCoeff = params.sparseTransform*gather(Im);
        
    % sparseCoeff(:,:,:, params.nsize(4)+1:end) = 5*sparseCoeff(:,:,:, params.nsize(4)+1:end);
    
    [sparseCoeffNorm, totalNorm] = params.sparseTransform.coeffNorm(sparseCoeff);
    sparseCoeff = params.sparseTransform.divideByNorm(sparseCoeff, sparseCoeffNorm, params.pNorm, params.l1Smooth);    
    
    % give temporal more weights
    if ( params.temporalScalingFactor ~= 1 )
        % sparseCoeff(:,:,:, params.nsize(4)+1:end) = params.temporalScalingFactor*sparseCoeff(:,:,:, params.nsize(4)+1:end);
        sparseCoeff = params.sparseTransform.scaleTemporalCoeff(sparseCoeff, params.temporalScalingFactor);
    end
    
    Im = gpuArray(params.sparseTransform'*sparseCoeff);
else
    sparseCoeff = params.sparseTransform*Im;

    if ( params.temporalScalingFactor ~= 1 )
        % sparseCoeff(:,:,:, params.nsize(4)+1:end) = params.temporalScalingFactor*sparseCoeff(:,:,:, params.nsize(4)+1:end);
        sparseCoeff = params.sparseTransform.scaleTemporalCoeff(sparseCoeff, params.temporalScalingFactor);
    end
    
    [sparseCoeffNorm, totalNorm] = params.sparseTransform.coeffNorm(sparseCoeff);
    sparseCoeff = params.sparseTransform.divideByNorm(sparseCoeff, sparseCoeffNorm, params.pNorm, params.l1Smooth);
    Im = params.sparseTransform'*sparseCoeff;
end

grad = 1/params.fftScaling*fftshift(fft2(ifftshift(Im))); 

function grad = gTikhonov(yu,params)
% compute gradient of the tikhonov operator
grad = 2*yu;

function grad = gTV(yu,params)
% compute gradient of TV operator

p = params.pNorm;

tempY = yu;
Im = params.fftScaling*fftshift(ifft2(ifftshift(tempY)));
Dx = params.TV*Im;

G = p*Dx.*(Dx.*conj(Dx) + params.l1Smooth).^(p/2-1);
grad = params.TV'*G;
% grad = fft2c(grad);
grad = 1/params.fftScaling*fftshift(fft2(ifftshift(grad)));
