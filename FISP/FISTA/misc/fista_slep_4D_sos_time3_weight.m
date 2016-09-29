function [x, funVal,  ValueL, mse]=fista_slep_4D(A, AT, y, lambda, opts)
%
%%
% Function fista_slep_4D
%      FISTA for 4D sum of squares reconstruction with time
%
%% Problem
%
%  min  1/2 \|Ax -y\|_F^2 + lambda * sum_i \| x(:,:,i).* [w(:,i), w(:,i), ...w(:,i)] \|_{2,1}
%                         + 0.5*rho* sum_i \|x(:,:,i)-x(:,:,i+1)\|_F^2
%
%  x(:,:,i)=[x^i_1, x^i_2, ..., x^i_c]
%           x^i_j is a concatenated vector from the 3D tensor
%           c is the number of coils
%           i goes from 1 to t, with t denoting the number of phases
%  y(index_b:index_e,:)=[y^i_1, y^i_2, ..., y^i_c]
%           y^i_j is a concatenated vector in the undersampled k-space
%           index_b:index_e corresponds to the size of y^i_j
%  the defulat value for w are all ones
%  x is 3D
%  y is 2D
%  w is 2D
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
% Composed by Jun Liu on June 06, 2011

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
% it is 3D tensor

if(isfield(opts,'rFlag'))
    if (opts.rFlag~=0 && opts.rFlag~=1)
        error('opts.rFlag error');
    end
else
    opts.rFlag=0; % no lambda_max
end

if(~isfield(opts,'w'))
    w=ones(n,t);
else
    w=opts.w;
    if norm( size(w) - [n,t] ) ~=0
        error('\n opts.w is incorrect');
    end
end

if (opts.rFlag)
    for i=1:t
        %temp=sum(abs(ATy(:,:,i)),2);
        
        temp=sqrt(sum(abs(ATy(:,:,i)).^2,2) );
        
        lambda_max(i,1)=max(temp); 
        % the lambda_max is for the case w=ones(n,t)
    end
else
    lambda_max=ones(t,1);
end
% for each time point, we compute its only maximal value

lambda_max=lambda_max./mean(w,2)';

lambda=lambda*lambda_max;
% it is a vector (not a scalar)

if(isfield(opts,'x0'))
    x0=opts.x0;
else
    x0=zeros(n,c,t);
end

if(isfield(opts,'L'))
    L=opts.L;
else
    L=1; % We assume that the maximum eigenvalue of A'A is over 1
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


proximity_operator = @(x, w, th) (soft4D(x, w, th));


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
        
        L=1;
        
        eta=0.6;
        lambda = lambda0;
    end
else
    opts.quickWarm=0; % default, do not use this
end


x=x0;

% compute A x
Ax=A(x);
% a 2D tensor

% assign xp with x, and Axp with Ax
xp=x; Axp=Ax; xxp=zeros(n,c,t);

% alphap and alpha are used for computing the weight in forming search point
alphap=0; alpha=1;

for iterStep=1:opts.maxIter
    fprintf('\n fista_slep_4D_sos_time iteration number: %d', iterStep);
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
    
    if (t==2)
        g(:,:,1)=g(:,:,1)+ rho* ( s(:,:,1)-s(:,:,2) );
        g(:,:,2)=g(:,:,2)+ rho* ( s(:,:,2)-s(:,:,1) );
    else
        g(:,:,1)=g(:,:,1)+ rho* ( s(:,:,1)-s(:,:,2) );
        g(:,:,t)=g(:,:,t)+ rho* ( s(:,:,t)-s(:,:,t-1) );
    end
    
    for i=2:(t-1)
        g(:,:,i)=g(:,:,i)+ rho* ( 2 * s(:,:,i)-s(:,:,i-1) -s(:,:,i+1));
    end
%       g(:,:, 2: (t-1) ) = g(:,:, 2: (t-1) ) + ...
%           rho * ( 2 * s(:,:, 2: (t-1) ) - s(:,:,1: (t-2)) - s(:,:,3: t) );
    
    % copy x and Ax to xp and Axp
    xp=x;    Axp=Ax;
    
    while (1)
        % let s walk in a step in the antigradient of s to get v
        % and then do the l1-norm regularized projection
        v=s-g/L;
        
        x = proximity_operator(v, w, lambda/L);
        
        v=x-s;  % the difference between the new approximate solution x
        % and the search point s
        
        % compute A x
        Ax=A(x);
        
        Av=Ax -As;
        r_sum=norm(v(:))^2; l_sum=norm(Av(:))^2;
        
        if (t==2)
            l_sum=l_sum + rho* ( norm (x(:,:,1)-x(:,:,2), 'fro')^2 - norm (s(:,:,1)-s(:,:,2), 'fro')^2 );
        else            
            l_sum=l_sum + rho* ( norm (x(:,:,1)-x(:,:,2), 'fro')^2 - norm (s(:,:,1)-s(:,:,2), 'fro')^2 );
            l_sum=l_sum + rho* ( norm (x(:,:,t)-x(:,:,t-1), 'fro')^2 - norm (s(:,:,t)-s(:,:,t-1), 'fro')^2 );
        end
        
        for i=2:(t-1)
            l_sum=l_sum + rho* ( norm (x(:,:,i)-x(:,:,i-1), 'fro')^2 - norm (s(:,:,i)-s(:,:,i-1), 'fro')^2 );
        end
        
%         if (r_sum <=1e-20)
%             bFlag=1; % this shows that, the gradient step makes little improvement
%             break;
%         end
        
        % the condition is ||Av||_F^2 <= (L) * ||v||_F^2
        if(l_sum <= r_sum * L)
            break;
        else
            L=max(1.5*L, l_sum/r_sum);
            fprintf('\n L=%e',L);
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
            funVal(iterStep)=funVal(iterStep) +lambda(i) * sum( w(:,i).* sqrt( sum(abs(x(:,:,i)).^2,2) ) );
        end
        
        if (t==2)
            funVal(iterStep)=funVal(iterStep) +0.5 *rho* ( norm (x(:,:,1)-x(:,:,2), 'fro')^2 );
        else
            funVal(iterStep)=funVal(iterStep) +0.5 *rho* ( norm (x(:,:,1)-x(:,:,2), 'fro')^2 );      
            funVal(iterStep)=funVal(iterStep) +0.5 *rho* ( norm (x(:,:,t)-x(:,:,t-1), 'fro')^2 );
        end
        
        for i=2:(t-1)
            funVal(iterStep)=funVal(iterStep) +0.5 *rho* ( norm (x(:,:,i)-x(:,:,i-1), 'fro')^2 );
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
                        funVal(iterStep)=funVal(iterStep) +lambda(i) * sum( w(:,i).* sqrt( sum(abs(x(:,:,i)).^2,2) ) );
                    end
                    
                    if (t==2)
                        funVal(iterStep)=funVal(iterStep) +rho* ( norm (x(:,:,1)-x(:,:,2), 'fro')^2 );
                    else
                        funVal(iterStep)=funVal(iterStep) +rho* ( norm (x(:,:,1)-x(:,:,2), 'fro')^2 );
                        funVal(iterStep)=funVal(iterStep) +rho* ( norm (x(:,:,t)-x(:,:,t-1), 'fro')^2 );
                    end
                    
                    for i=2:(t-1)
                        funVal(iterStep)=funVal(iterStep) +rho* ( norm (x(:,:,i)-x(:,:,i-1), 'fro')^2 );
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
                        funVal(iterStep)=funVal(iterStep) +lambda(i) * sum( w(:,i).* sqrt( sum(abs(x(:,:,i)).^2,2) ) );
                    end
                    
                    if (t==2)
                        funVal(iterStep)=funVal(iterStep) +rho* ( norm (x(:,:,1)-x(:,:,2), 'fro')^2 );
                    else
                        funVal(iterStep)=funVal(iterStep) +rho* ( norm (x(:,:,1)-x(:,:,2), 'fro')^2 );
                        funVal(iterStep)=funVal(iterStep) +rho* ( norm (x(:,:,t)-x(:,:,t-1), 'fro')^2 );
                    end
                    
                    for i=2:(t-1)
                        funVal(iterStep)=funVal(iterStep) +rho* ( norm (x(:,:,i)-x(:,:,i-1), 'fro')^2 );
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

