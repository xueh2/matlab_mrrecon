function [x, funVal,  ValueL, mse]=fista_slep_1D_wavelet(A, AT, y, lambda, opts)
%
%%
% Function fista_slep_1D_wavelet
%      FISTA for 1D reconstruction
%
%% Problem
%
%  min  1/2 \|Ax -y\|_F^2 + lambda * sum_i \| W x(:,i)\|_{1}
%
%  x=[x_1, x_2, ..., x_t]
%           x_i is a 1D vector
%           i goes from 1 to t, with t denoting the number of points
%  y=[y_1; y_2; ...; y_t]
%           y_i is a 1D vector, y_i=[y_i^1, y_i^2, ..., y_i^c]
%  x is 1D
%  y is 1D
%
%
%  Note: \|x\|_1 is defined as the summation of the magnitude of each entry
%  (which is complex)
%
%% Input parameters:
%
%  A, AT-      matrix
%  y -         response data
%
%% Output parameters:
%  x-         Solution
%  funVal-    Function value during iterations
%

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

t=numel(y);
% the number of time points

if (t<2)
    error('\n t should be over 2');
end

% ATy
ATy=AT*y;

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
    x0=zeros(n,1);
end

if(isfield(opts,'L'))
    L=opts.L;
else
    L=0.1; % we assume that the maximum eigenvalue of A'A is over 1
end

if(~isfield(opts,'maxIter'))
    opts.maxIter=1000;
end

if(~isfield(opts,'wavelet'))
    opts.wavelet=1;
end

if(~isfield(opts,'tFlag'))
    opts.tFlag=3;
end

if(~isfield(opts,'tol'))
    opts.tol=1e-5;
end

keepLL = 0;
if(isfield(opts,'keepLL'))
    keepLL=opts.keepLL;
else
    keepLL = 0;
end

W=opts.W;
WT=opts.WT;
scale=opts.scale;

p=zeros(length(W(x0)),1);
q=p;
proximity_operator = @(x, th) Dykstra(x,th);


WATy= W(ATy);

if (opts.rFlag)
    lambda_max=max(abs(WATy(:))); % we use a unified lambda_max
    lambda=lambda*lambda_max;
else
    lambda_max=opts.lambda_max;
    lambda=lambda_max;
end
% for each time point, we compute its only maximal value

if (opts.verbose)
    lambda_max
end
% rho=rho* lambda_max;
% it is a scalar

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
    else
        lambda_bar=lambda;
    end
else
    lambda_bar=lambda;
    
    opts.quickWarm=0; % default, do not use this
end


x=x0;

% compute A x
Ax=A*x;
% a 2D matrix

% assign xp with x, and Axp with Ax
xp=x; Axp=Ax; xxp=zeros(n,1);

% alphap and alpha are used for computing the weight in forming search point
alphap=0; alpha=1;

for iterStep=1:opts.maxIter
    
    if (opts.verbose)
        fprintf('\n fista_slep_2D_csm iteration number: %d', iterStep);
    end
    % --------------------------- step 1 ---------------------------
    % compute search point s based on xp and x (with beta)
    beta=(alphap-1)/alpha;    s=x + beta* xxp;
    
    % --------------------------- step 2 ---------------------------
    % line search for L and compute the new approximate solution x
    
    % compute the gradient (g) at s
    As=Ax + beta* (Ax-Axp);
    
    % compute AT As
    ATAs=AT*As;
    
    % obtain the gradient g
    g=ATAs-ATy;
    g=g(:);
        
    % copy x and Ax to xp and Axp
    xp=x;    Axp=Ax;
    
    while (1)
        % let s walk in a step in the antigradient of s to get v
        % and then do the l1-norm regularized projection
        v=s-g/L;
        
        x = proximity_operator(v, lambda/L);
        
        v=x-s;  % the difference between the new approximate solution x
        % and the search point s
        
        % compute A x
        Ax=A*x;
        
        Av=Ax -As;
        r_sum=norm(v(:))^2; l_sum=norm(Av(:))^2;
        
%         if (r_sum <=1e-20)
%             bFlag=1; % this shows that, the gradient step makes little improvement
%             break;
%         end

        %if (r_sum<=1e-20)
        if (r_sum<=1e-5)
            break;
        end
        
        % the condition is ||Av||_F^2 <= (L) * ||v||_F^2
        if(l_sum <= r_sum * L)
            break;
        else
            L=max(1.5*L, l_sum/r_sum);
            if (opts.verbose)
                fprintf('\n L=%e, l_sum=%e, r_sum * L=%e, diff=%e',L, l_sum, r_sum * L, l_sum - r_sum * L);
            end
        end
        
    end
    
    if(nargout>=3)
        ValueL(iterStep)=L;
    end
    
    % --------------------------- step 3 ---------------------------
    % update alpha and alphap, and check whether converge
    alphap=alpha; alpha= (1+ sqrt(4*alpha*alpha +1))/2;
    
    xxp=x-xp;   Axy=Ax-y;
    
    error_data_fidelity=norm(Axy(:))^2/2;
    if(nargout>=2)
        funVal(iterStep)= error_data_fidelity+ lambda* sum(abs( W(x) ) ); 
    end
    
    if (opts.verbose)
        fprintf('\n x: %e, %e', max( abs(x) ), sum( abs(x) ));
        fprintf('   %e, %e, %e',error_data_fidelity, funVal(iterStep)- error_data_fidelity, funVal(iterStep));
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
                    funVal(iterStep)=norm(Axy(:))^2/2 + lambda* sum(abs( W(x) ) ); 

                    %fprintf('iterStep: %d, lambda =%e, function value = %e', iterStep, lambda, funVal(iterStep));
                end
                
                if (abs( funVal(iterStep) - funVal(iterStep-1) ) <= opts.tol)
                    break;
                end
            end
        case 1
            if iterStep>=2
                
                if(nargout<2)
                    funVal(iterStep)=norm(Axy(:))^2/2 + lambda* sum(abs( W(x) ) ); 
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
            
            if iterStep>=20
                norm_xxp=norm(xxp(:),'fro');
                if ( norm_xxp <=opts.tol)
                    break;
                end
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


function x=Dykstra(v,th)
    
    % v: x-t domain
    
    
    Wv=W(v);
    % wavelet coefficient
    
    
    if (lambda > lambda_bar )
        p=0;
        q=0;
        x=Wv;
    else
        p=sign(p).*min(abs(p),th);
        x=Wv-p-q;
    end
    
    for i=1:10
        if (opts.wavelet==1)
            xplusp=x+p;
            if (keepLL)
                yy=soft_keep_LL(xplusp,th); 
            else
                yy=soft(xplusp,th); 
            end
            yy=yy(:);
            %yy=soft_temporal_3d_t(xplusp,th, scale); % the difference between the first and last temporal phase is not penalalized
        else
            if (keepLL)
                yy=soft_keep_LL(x+p,th);
            else
                yy=soft(x+p,th);
            end
        end
        
        p=p+(x-yy);
        
        x=W( WT(yy+q) );
        q=q+(yy-x);
                
        %fprintf('\n %e', 0.5* norm(x-Wv).^2 + th * sum(abs(x)));        
        
        if norm(x-yy)< norm(x)* 1e-4
            break;
        end
        
    end
    
    x=WT(x);
end
end