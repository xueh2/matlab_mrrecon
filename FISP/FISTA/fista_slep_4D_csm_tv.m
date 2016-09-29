function [x, funVal,  ValueL, mse]=fista_slep_4D_csm_tv(A, AT, y, lambda, opts)
%
%%
% Function fista_slep_4D_csm_tv
%      FISTA for 4D MRA reconstruction with csm and tv for temporal
%
%% Problem
%
%  min  1/2 \|Ax -y\|_F^2 + lambda * sum_i \| x(:,i)\|_{1}
%                         + 0.5*rho* sum_i \|x(:,i)-x(:,i+1)\|_{1}
%
%  x=[x_1, x_2, ..., x_t]
%           x_i is a concatenated vector from the 3D tensor
%           i goes from 1 to t, with t denoting the number of phases
%  y=[y_1; y_2; ...; y_t]
%           y_i is a concatenated vector in the undersampled k-space of
%           phase i, y_i=[y_i^1, y_i^2, ..., y_i^c]
%  x is 2D
%  y is 2D
%
%  Note: here x is complex, and \|x\|_1 is defined as
%                               \|real(x)\|_1 + \|imag(x)\|_1
%        This is different from the one used in fista_slep_4D_csm
%
%
%% Input parameters:
%
%  A, AT-      partial fourier operator
%  y -         k-space data
%
%% Output parameters:
%  x-         Solution
%  funVal-    Function value during iterations
%
%%
% Composed by Jun Liu on June 15, 2011

%% Verify and initialize the parameters
%%

% Verify the number of input parameters
if (nargin <5)
    error('\n Inputs: A, AT, y, lambda, opts should be specified!\n');
end

if(isfield(opts,'n'))
    n=opts.n;
else
    error('opts.n should be specified');
end
% n denotes the number of rows of x

if(isfield(opts,'rho'))
    rho=opts.rho;
else
    rho=0;
end
% rho is the parameter for fusing different time...

global sizeImage;
% it contains four elements [1d,2d,3d,c, t]

c=sizeImage(4);
% the number of coils
t=sizeImage(5);
% the number of time phase

if (t<2)
    error('\n t should be over 2');
end

% ATy
ATy=AT(y);
% it is 2D matrix

if(isfield(opts,'rFlag'))
    if (opts.rFlag~=0 && opts.rFlag~=1)
        error('opts.rFlag error');
    end
else
    opts.rFlag=0; % no lambda_max
end

if(isfield(opts,'x0'))
    x0=opts.x0;
else
    x0=zeros(n,t);
end

if(isfield(opts,'L'))
    L=opts.L;
else
    L=0.1; % we assume that the maximum eigenvalue of A'A is over 1
end

if(~isfield(opts,'maxIter'))
    opts.maxIter=1000;
end

if(~isfield(opts,'tFlag'))
    opts.tFlag=3;
end

if(~isfield(opts,'tol'))
    opts.tol=1e-5;
end


if (opts.rFlag)
    lambda_max= max( max(real(ATy(:))), max(imag(ATy(:))) )
    % we use a unified lambda_max
    % for all temporal phases, and the weights for real and imaginary part
    % are same. 
else
    lambda_max=1;
end
% for each time point, we compute its only maximal value

rho=rho* lambda_max;
lambda=lambda*lambda_max;
% it is a vector (not a scalar)

bFlag=0; % this flag tests whether the gradient step only changes a little

if (isfield(opts,'quickWarmStart'))
    if (opts.quickWarmStart)
        
        if (~isfield(opts,'rFlag'))
            error('\n opts.rFlag should equal to 1');
        else
            if (opts.rFlag~=1)
                error('\n opts.rFlag should equal to 1');
            end
        end
                
        lambda0 = 0.5*lambda_max;
        
        if(isfield(opts,'lambda_bar'))
            lambda_bar=opts.lambda_bar*lambda_max;
        else            
            lambda_bar = 0.001*lambda_max;
        end
        
        L=0.1;
        
        eta=0.6;
        lambda = lambda0;
    end
else
    opts.quickWarm=0; % default, do not use this
end


x=x0;

% compute A x
Ax=A(x);
% a 2D matrix

% assign xp with x, and Axp with Ax
xp=x; Axp=Ax; xxp=zeros(n,t);

% initialize z0 for real and imaginary parts
z0_real=zeros(n,t-1);
z0_imag=zeros(n,t-1);

% alphap and alpha are used for computing the weight in forming search point
alphap=0; alpha=1;

for iterStep=1:opts.maxIter
    fprintf('\n fista_slep_4D_csm iteration number: %d', iterStep);
    % --------------------------- step 1 ---------------------------
    % compute search point s based on xp and x (with beta)
    beta=(alphap-1)/alpha;    s=x + beta* xxp;
    
    % --------------------------- step 2 ---------------------------
    % line search for L and compute the new approximate solution x
    
    % compute the gradient (g) at s
    As=Ax + beta* (Ax-Axp);
    
    % compute AT As
    ATAs=AT(As);
    
    % obtain the gradient g
    g=ATAs-ATy;
    
    % copy x and Ax to xp and Axp
    xp=x;    Axp=Ax;
    
    while (1)
        % let s walk in a step in the antigradient of s to get v
        % and then do the l1-norm regularized projection
        v=s-g/L;
        v_real=real(v); v_imag=imag(v);
        
%         x_real = proximity_operator(v_real, lambda/L);
%         x_imag = proximity_operator(v_imag, lambda/L);
        
        [x_real, z_real, gap ]=tesla_proj(v_real,z0_real,lambda/L,rho/L,n,t,100,1e-10,1,6);
        [x_imag, z_imag, gap ]=tesla_proj(v_imag,z0_imag,lambda/L,rho/L,n,t,100,1e-10,1,6);
        z0_real=z_real;
        z0_imag=z_imag;
        x=complex(x_real,x_imag);
        
        v=x-s;  % the difference between the new approximate solution x
        % and the search point s
        
        % compute A x
        Ax=A(x);
        
        Av=Ax -As;
        r_sum=norm(v(:))^2; l_sum=norm(Av(:))^2;
        
%         if (r_sum <=1e-20)
%             bFlag=1; % this shows that, the gradient step makes little improvement
%             break;
%         end

        if (r_sum<=1e-20)
            break;
        end
        
        % the condition is ||Av||_F^2 <= (L) * ||v||_F^2
        if(l_sum <= r_sum * L)
            break;
        else
            L=max(1.5*L, l_sum/r_sum);
            fprintf('\n L=%e, l_sum=%e, r_sum * L=%e, diff=%e',L, l_sum, r_sum * L, l_sum - r_sum * L);
        end
    end
    
    if(nargout>=3)
        ValueL(iterStep)=L;
    end
    
    % --------------------------- step 3 ---------------------------
    % update alpha and alphap, and check whether converge
    alphap=alpha; alpha= (1+ sqrt(4*alpha*alpha +1))/2;
    
    xxp=x-xp;   Axy=Ax-y;
    
    if(nargout>=2)
        funVal(iterStep)=norm(Axy(:))^2/2;
        
        for i=1:t
            funVal(iterStep)=funVal(iterStep) + lambda * ( sum( real(x(:,i)) ) + sum( real(x(:,i)) ) );
        end
       
        for i=2:t
            funVal(iterStep)=funVal(iterStep) + rho* ( sum( real(x(:,i)-x(:,i-1) ) ) + sum( imag(x(:,i)-x(:,i-1) ) ) );
        end
    end
    
    if (opts.quickWarmStart)        
        lambda=max(eta*lambda,lambda_bar);
    end
    
    if(nargout>=4)
        if(~isfield(opts,'groundTruth'))
            error('groundTruth should be given');
        else
            groundTruth=opts.groundTruth;
            mse(iterStep)=norm(real(x)-groundTruth,2)/numel(groundTruth);
        end        
    end
    
    if (bFlag)
        % fprintf('\n The program terminates as the gradient step changes the solution very small.');
        break;
    end
    
    switch(opts.tFlag)
        case 0
            if iterStep>=2
                
                if(nargout<2)
                    funVal(iterStep)=norm(Axy(:),'fro')^2/2
                    
                    
                    for i=1:t
                        funVal(iterStep)=funVal(iterStep) + lambda * ( sum( real(x(:,i)) ) + sum( real(x(:,i)) ) );
                    end
                    
                    for i=2:t
                        funVal(iterStep)=funVal(iterStep) + rho* ( sum( real(x(:,i)-x(:,i-1) ) ) + sum( imag(x(:,i)-x(:,i-1) ) ) );
                    end

                    %fprintf('iterStep: %d, lambda =%e, function value = %e', iterStep, lambda, funVal(iterStep));
                end
                
                if (abs( funVal(iterStep) - funVal(iterStep-1) ) <= opts.tol)
                    break;
                end
            end
        case 1
            if iterStep>=2
                
                if(nargout<2)
                    funVal(iterStep)=norm(Axy(:),'fro')^2/2
                    
                    for i=1:t
                        funVal(iterStep)=funVal(iterStep) + lambda * ( sum( real(x(:,i)) ) + sum( real(x(:,i)) ) );
                    end
                    
                    for i=2:t
                        funVal(iterStep)=funVal(iterStep) + rho* ( sum( real(x(:,i)-x(:,i-1) ) ) + sum( imag(x(:,i)-x(:,i-1) ) ) );
                    end
                    %fprintf('iterStep: %d, lambda =%e, function value = %e', iterStep, lambda, funVal(iterStep));
                end
                
                if (abs( funVal(iterStep) - funVal(iterStep-1) ) <=...
                        opts.tol* funVal(iterStep-1))
                    break;
                end
            end
        case 2
            if ( funVal(iterStep)<= opts.tol)
                break;
            end
        case 3
            norm_xxp=norm(xxp(:),'fro');
            if ( norm_xxp <=opts.tol)
                break;
            end
        case 4
            norm_xp=norm(xp(:),'fro');    norm_xxp=norm(xxp(:),'fro');
            if ( norm_xxp <=opts.tol * max(norm_xp,1))
                break;
            end
        case 5
            if iterStep>=opts.maxIter
                break;
            end
    end
end

