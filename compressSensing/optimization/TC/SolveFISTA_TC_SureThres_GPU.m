% Copyright ©2010. The Regents of the University of California (Regents). 
% All Rights Reserved. Contact The Office of Technology Licensing, 
% UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620, 
% (510) 643-7201, for commercial licensing opportunities.

% Authors: Arvind Ganesh, Allen Y. Yang, Zihan Zhou.
% Contact: Allen Y. Yang, Department of EECS, University of California,
% Berkeley. <yang@eecs.berkeley.edu>

% IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, 
% SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, 
% ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF 
% REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, 
% PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO 
% PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

%% This function is modified from Matlab code proximal_gradient_bp

function [x_hat,nIter, costValue] = SolveFISTA_TC_SureThres_GPU(A,AT,b,W,WT, varargin)

% b - m x 1 vector of observations/data (required input)
% A - m x n measurement matrix (required input)
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
% maxIter - maxilambdam number of iterations
%         - DEFAULT 10000, if omitted or -1.
% lineSearchFlag - 1 if line search is to be done every iteration
%                - DEFAULT 0, if omitted or -1.
% continuationFlag - 1 if a continuation is to be done on the parameter lambda
%                  - DEFAULT 1, if omitted or -1.
% eta - line search parameter, should be in (0,1)
%     - ignored if lineSearchFlag is 0.
%     - DEFAULT 0.9, if omitted or -1.
% lambda - relaxation parameter
%    - ignored if continuationFlag is 1.
%    - DEFAULT 1e-3, if omitted or -1.
% outputFileName - Details of each iteration are dumped here, if provided.
%
% x_hat - estimate of coeeficient vector
% numIter - number of iterations until convergence
%
%
% References
% "Robust PCA: Exact Recovery of Corrupted Low-Rank Matrices via Convex Optimization", J. Wright et al., preprint 2009.
% "An Accelerated Proximal Gradient Algorithm for Nuclear Norm Regularized Least Squares problems", K.-C. Toh and S. Yun, preprint 2009.
%
% Arvind Ganesh, Summer 2009. Questions? abalasu2@illinois.edu

startT = tic;

sizeInfo = size(b);

x0 = parallel.gpu.GPUArray.zeros(sizeInfo) ;
xG = [];

STOPPING_GROUND_TRUTH = -1;
STOPPING_DUALITY_GAP = 1;
STOPPING_SPARSE_SUPPORT = 2;
STOPPING_OBJECTIVE_VALUE = 3;
STOPPING_SUBGRADIENT = 4;
STOPPING_DEFAULT = STOPPING_SUBGRADIENT;

% Default values for optional arguments.
rho = 1;
stoppingCriterion = STOPPING_DEFAULT;
isInit = true;
ratioFISTA = 1.5;

% Parse the optional inputs.
if (mod(length(varargin{:}), 2) ~= 0 ),
    error(['Extra Parameters passed to the function ''' mfilename ''' lambdast be passed in pairs.']);
end
parameterCount = length(varargin{:})/2;

varargin = varargin{1};

for parameterIndex = 1:parameterCount,
    parameterName = varargin{parameterIndex*2 - 1};
    parameterValue = varargin{parameterIndex*2};
    switch lower(parameterName)
        case 'stoppingcriterion'
            stoppingCriterion = parameterValue;
        case 'groundtruth'
            xG = parameterValue;
        case 'tolerance'
            tolerance = parameterValue;
        case 'linesearchflag'
            lineSearchFlag = parameterValue;
        case 'lambda'
            lambda_bar_factor = parameterValue;
        case 'maxiteration'
            maxIter = parameterValue;
        case 'isnonnegative'
            isNonnegative = parameterValue;
        case 'continuationflag'
            continuationFlag = parameterValue;
        case 'initialization'
            xk = parameterValue;
            if ~all(size(xk)==sizeInfo)
                error('The dimension of the initial xk does not match.');
            end
        case 'eta'
            eta = parameterValue;
            if ( eta <= 0 || eta >= 1 )
                disp('Line search parameter out of bounds, switching to default 0.9') ;
                eta = 0.9 ;
            end
        case 'csm'
            coilSensitivity = parameterValue;
            csm_sos = sum(abs(coilSensitivity).^2, 3);
            rho = max(abs(csm_sos(:)));
        case 'tv'
            useTV = parameterValue;
        case 'startparam'
            isInit = false;
            t_km1 = parameterValue.t_km1;
            t_k = parameterValue.t_k;
            xkm1 = parameterValue.xkm1;
            L = parameterValue.L;
            lambda = parameterValue.lambda;
            lambda_bar = parameterValue.lambda_bar;
        case 'ratiofista'
            ratioFISTA = parameterValue;
        case 'gfactor'
            gFactor = parameterValue;
        case 'mocoer'
            mocoer = parameterValue;
        case 'senmap'
            senMap = parameterValue;
            conjSenMap = conj(senMap);
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end
clear varargin

if ( ~isvar('gFactor') )
    gFactor = [];
end

if ( ~isvar('mocoer') )
    mocoer = [];
end

if ( ~isvar('senMap') )
    senMap = [];
end

if stoppingCriterion==STOPPING_GROUND_TRUTH && isempty(xG)
    error('The stopping criterion must provide the ground truth value of x.');
end

if ~isvar('maxIter')
    maxIter = 200 ;
end
if ~isvar('tolerance');
    tolerance = 1e-6;
end
if ~isvar('lambda_bar_factor');
    lambda_bar_factor = 0.01;
end

G = @(x) AT(A(x)); %function handle

%% Initializin optimization variables

aSize = size(b);
sx = aSize(1);
sy = aSize(2);
fftScaling = sqrt(sx*sy);

if ~isvar('xk')
    xk = zeros(sizeInfo) ;
end
c = AT(b) ; % function handle
% cc = ifft2c(c);
cc = fftScaling*fftshift(ifft2(ifftshift(c)));
nIter = 0 ;
keep_going = 1 ;
eta = 1;
beta = 2;

if isInit
    t_k = 1 ; 
    t_km1 = 1 ;
    L0 = 0.1 ;
    L = L0 / rho;

    % norminf = mean(abs(cc(:)))/ratioFISTA;
    
    cc = SoS_Image_TemporalArray(gather(cc));
    cc = W*double(cc+i*cc);    
    norminf = max(abs(real(cc(:))))*lambda_bar_factor;
    
    lambda0 = norminf;
    lambda_bar = lambda_bar_factor*L*norminf ;
    lambda = lambda0 ;

    nz_x = (abs(xk)> eps*10);
    xkm1 = xk;
end

while keep_going && (nIter < maxIter)

    nIter = nIter + 1 ;
    yk = xk + ((t_km1-1)/t_k)*(xk-xkm1) ;    
    stop_backtrack = 0 ;
   
    gradient = G(yk) - c; % gradient of f at yk

    while ~stop_backtrack
        
        gk = yk - (1/(L * rho))*gradient ;
        % xkp1 = proximity_operator(gk, lambda/(L * rho));

        % compute the sure threshold for gk
        gkIm = double(gather( fftScaling*fftshift(ifft2(ifftshift(gk))) ));
        gkIm = SoS_Image_TemporalArray(gkIm);
        gkCoeff = W*(gkIm+i*gkIm);
        
        if ( numel(size(gkCoeff) ) == 4 )
            finestCoeff = real(gkCoeff(:,:,end));
        else
            finestCoeff = real(gkCoeff(sizeInfo(1)+1:end, sizeInfo(2)+1:end, :));
        end
        
        noiseSigma = median(abs(finestCoeff(:)))*1.4286;
        thres = EstimateSureShrinkThreshold(abs(real(gkCoeff)), noiseSigma);
        xkp1 = proximity_operator(gk, 0.5*(lambda/(L * rho)+thres/ratioFISTA));
        % xkp1 = proximity_operator(gk, lambda/(L * rho));
        
        temp1 = funF(xkp1); % f(xkp1)
        % temp2 = funF(yk) + gradient(:)'*(xkp1(:)-yk(:)) + ((L * rho)/2)*norm(xkp1(:)-yk(:),'fro')^2 ; % f(yk) + (xkp1 - yk)'*nabla(f)(yk) + (L/2)*norm(xkp1-yk)^2
        v1 = xkp1-yk;
        temp2 = funF(yk) + gradient(:)'*(xkp1(:)-yk(:)) + ((L * rho)/2)*abs(v1(:)'*v1(:)) ;
        
        if abs(temp1) <= abs(temp2)
            stop_backtrack = 1 ;
        else
            L = L*beta ;
        end
    end

    switch stoppingCriterion
        case STOPPING_GROUND_TRUTH
            keep_going = norm(xG(:)-xkp1(:))>tolerance;
        case STOPPING_SUBGRADIENT
            sk = L*rho*(yk-xkp1) + G(xkp1-yk);
            % FIXME: norm or norm('fro') for complexes?
            % lhs_sg = norm(sk(:));
            lhs_sg = normGPU(sk);
            % rhs_sg = L*rho*max(1,norm(xkp1(:)));
            rhs_sg = L*rho*max(1,normGPU(xkp1));
            keep_going = ((lhs_sg / rhs_sg) > tolerance) ;
            disp(['stopcriterion(' num2str(nIter) ') ' num2str(abs(gather(lhs_sg / rhs_sg)), 15) ' < ' num2str(tolerance, 15)]);
        case STOPPING_SPARSE_SUPPORT
            % compute the stopping criterion based on the change
            % of the number of non-zero components of the estimate
            nz_x_prev = nz_x;
            nz_x = (abs(xkp1)>eps*10);
            num_nz_x = sum(nz_x(:));
            num_changes_active = (sum(nz_x(:)~=nz_x_prev(:)));
            if num_nz_x >= 1
                criterionActiveSet = num_changes_active / num_nz_x;
                keep_going = (criterionActiveSet > tolerance);
            end
        case STOPPING_OBJECTIVE_VALUE
        case STOPPING_DUALITY_GAP
            error('Duality gap is not a valid stopping criterion for PGBP.');
        otherwise
            error('Undefined stopping criterion.');
    end
    
    lambda = max(eta*lambda,lambda_bar) ;
    t_kp1 = 0.5*(1+sqrt(1+4*t_k*t_k)) ;
    
    t_km1 = t_k ;
    t_k = t_kp1 ;
    xkm1 = xk ;
    xk = xkp1;
    
    % L = L0;
end
% Wx = W*double(gather(ifft2c(xk)));
Wx = W*double(gather(fftScaling*fftshift(ifft2(ifftshift(xk)))));
residual=b-A(xk);
% costValue = 0.5*norm(residual(:))^2 + lambda0 * norm(Wx(:),1);
costValue = 0.5*residual(:)'*residual(:) + lambda0 * sum(abs(Wx(:)));

x_hat = xk;

x_hat = gather(x_hat);
costValue = gather(costValue);

disp(['Time cost is ' num2str(toc(startT))]);

    %% norm computation function
    function v = funF(x)
        y = A(x) - b;
        v = 0.5 * abs(y(:)'*y(:)); 
    end

    function v = normGPU(x)
        v = sqrt(abs(x(:)'*x(:))); 
    end

    function v = proximity_operator(x, th)
        x2 = double(gather( fftScaling*fftshift(ifft2(ifftshift(x))) ));
        
        if ( ~isempty(senMap) )
            x2 = sum(conjSenMap.*x2, 3);
        end
        
        if ( ~isempty(mocoer) )
            x2 = mocoer*(1e5.*x2);
            x2 = x2 ./ 1e5;
        end
        
        x2 = W*x2;
%         if ( isvarname('gFactor') )
%             x2 = W.softThresh(x2,gather(th),gFactor);
%         else
%             x2 = W.softThresh(x2,gather(th),[]);
%         end
        x2 = W.softThresh(x2,gather(th));
        % x2 = soft(x2,double(gather(th)));
        x2 = WT*x2;
        
        if ( ~isempty(mocoer) )
            x2 = mocoer'*(1e5.*x2);
            x2 = x2 ./ 1e5;
        end
        
        if ( ~isempty(senMap) )
            x2 = repmat(x2, [1 1 size(x, 3) 1]);
            x2 = x2 .* senMap;
        end
        
        % v = fft2c(gpuArray(x2));
        v = 1/fftScaling*fftshift(fft2(ifftshift(gpuArray(x2))));
    end
end