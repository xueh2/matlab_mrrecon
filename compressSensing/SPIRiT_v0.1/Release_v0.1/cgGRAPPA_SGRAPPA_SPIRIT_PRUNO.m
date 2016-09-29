function [res, RESVEC] = cgGRAPPA_SGRAPPA_SPIRIT_PRUNO(y, ygrappa, SGRA, GOP, PRU, weighting, nIter, thres, x0)
%
%
%  res = cgGRAPPA_SGRAPPA_SPIRIT_PRUNO(y, ygrappa, SGRA, GOP, PRU, nIter, lambda, x0)
%  
%  Implementation of the Cartesian, union kernel reconstruction
%  Solve grappa, self-consistent grappa, spirit, pruno union recon
%
%  Input:
%		y -	Undersampled k-space data. Make sure that empty entries are zero
%			or else they will not be filled.
%       ygrappa - grappa reconed kspace
%       SGRA - self-consistent grappa operator
%       GOP - SPIRiT operator
%		PRU -	the PURNO operator obtained by calibration, implementing N'N
%		nIter -	Maximum number of iterations
%		lambda-	Tykhonov regularization parameter
%       x0    - Initil value
%
% Outputs:
%		res - Full k-space matrix
%
%
% Example:
%
%   Hui Xue, July 2011
%

t=tic;

if nargin < 4
	lambda = 0;
end

if nargin < 5
     x0 = y;
end

[sx,sy,nCoils] = size(y);

idx_acq = find(abs(y)>0);
idx_nacq = find(abs(y)==0);

% solve A*x = b
% compute - b
b = [];
bGrappa = [];
bSGrappa = [];
bSPIRiT = [];
bPRUNO = [];

% grappa
if ( ~isempty(ygrappa) )
    bGrappa = ygrappa - y;
    bGrappa = bGrappa(:);
    b = weighting(1) * bGrappa;
end

% s-grappa
if ( ~isempty(SGRA) )
    bSGrappa = y;
    bSGrappa = bSGrappa(:);
    b = [b; weighting(2) * bSGrappa];
end

% spirit
if ( ~isempty(GOP) )
    bSPIRiT = GOP*y;
    bSPIRiT = -bSPIRiT(:);
    b = [b; weighting(3) * bSPIRiT];
end

% pruno
if ( ~isempty(PRU) )
    bPRUNO = PRU*y; % (N'N)*y 
    
    bPRUNO = bPRUNO(idx_nacq); % apply Dc
    
    bPRUNO = -bPRUNO(:);
    b = [b; weighting(4) * bPRUNO];
end

N = length(b);
if ( N == 0 )
    disp('No kernels are unioned ... ');
    res = y;
    RESVEC = 0;
    return;
end

[tmpres,FLAG,RELRES,ITER,RESVEC] = lsqr(@aprod,b,thres,nIter, [],[],x0(idx_nacq),ygrappa, SGRA, GOP, PRU,sx,sy,nCoils,idx_acq,idx_nacq,weighting);
FLAG
RELRES
ITER

res = y;
res(idx_nacq) = tmpres;

disp(['cgGRAPPA_SGRAPPA_SPIRIT_PRUNO : ' num2str(toc(t))]);

function [res,tflag] = aprod(x,ygrappa, SGRA, GOP, PRU,sx,sy,nCoils,idx_acq,idx_nacq,weighting,tflag)

	if strcmp(tflag,'transp');
        
        aN = sx*sy*nCoils;
        N = length(x);
        
        res = zeros(length(idx_nacq), 1);
        rGrappa = [];
        rSGrappa = [];
        rSPIRiT = [];
        rPRUNO = [];
        kernelUnioned = 0;
       
        % grappa
        if ( ~isempty(ygrappa) )
            tmpy = reshape(x(kernelUnioned*aN+1:(kernelUnioned+1)*aN),sx,sy,nCoils);
            rGrappa = tmpy(idx_nacq);
            res = res + weighting(1)* rGrappa(:);
            kernelUnioned = kernelUnioned + 1;
        end

        % s-grappa
        if ( ~isempty(SGRA) )
            tmpy = reshape(x(kernelUnioned*aN+1:(kernelUnioned+1)*aN),sx,sy,nCoils);
            rSGrappa = SGRA'*tmpy;
            rSGrappa = rSGrappa - tmpy;
            rSGrappa = rSGrappa(idx_nacq);
            res = res + weighting(2)* rSGrappa(:);
            kernelUnioned = kernelUnioned + 1;
        end

        % spirit
        if ( ~isempty(GOP) )
            tmpy = reshape(x(kernelUnioned*aN+1:(kernelUnioned+1)*aN),sx,sy,nCoils);
            rSPIRiT = GOP'*tmpy;
            rSPIRiT = rSPIRiT(idx_nacq);
            res = res + weighting(3)* rSPIRiT(:);
            kernelUnioned = kernelUnioned + 1;
        end

        % pruno
        if ( ~isempty(PRU) )
%             tmpy = reshape(x(kernelUnioned*aN+1:(kernelUnioned+1)*aN),sx,sy,nCoils);                       
%             rPRUNO = PRU*tmpy;
%             rPRUNO = rPRUNO(idx_nacq);
%             res = [res; rPRUNO(:)];
%             kernelUnioned = kernelUnioned + 1;
            
            tmpy = zeros(sx,sy,nCoils);
            tmpy(idx_nacq) = x(kernelUnioned*aN+1:end);
            rPRUNO = PRU*tmpy;
            rPRUNO(idx_acq) = 0;
            rPRUNO = rPRUNO(idx_nacq);
            res = res + weighting(4)* rPRUNO(:);
        end
	
    else
        res = [];
        rGrappa = [];
        rSGrappa = [];
        rSPIRiT = [];
        rPRUNO = [];

		tmpx = zeros(sx,sy,nCoils);
		tmpx(idx_nacq) = x;
        
        % grappa
        if ( ~isempty(ygrappa) )
            rGrappa = tmpx;
            res = [res; weighting(1)* rGrappa(:)];
        end

        % s-grappa
        if ( ~isempty(SGRA) )
            rSGrappa = SGRA*tmpx;
            rSGrappa = rSGrappa - tmpx;
            res = [res; weighting(2)* rSGrappa(:)];
        end

        % spirit
        if ( ~isempty(GOP) )
            rSPIRiT = GOP*tmpx;
            res = [res; weighting(3)* rSPIRiT(:)];
        end

        % pruno
        if ( ~isempty(PRU) )
            rPRUNO = PRU*tmpx;            
            rPRUNO = rPRUNO(idx_nacq);            
            res = [res; weighting(4)* rPRUNO(:)];            
        end
	end
