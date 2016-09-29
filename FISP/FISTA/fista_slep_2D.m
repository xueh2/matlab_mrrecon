function [x, funVal,  ValueL, mse]=fista_slep_2D(A, AT, y, lambda, opts)
%
%%
% Function fista_slep_2D
%      FISTA for 2D sum of squares reconstruction
%
%% Problem
%
%  min  1/2 \|Ax -y\|_F^2 + lambda *  \| x.* [w, w, ..., w] \|_{2,1}
%
%  x=[x(:,1), x(:,2), ..., x(:,c)]
%           x(:,i) is a concatenated vector from the 2D matrix
%           c is the number of coils
%  y=[y_1; y_2; y_c]
%           y_i is a concatenated vector in the undersampled k-space
%  the defulat value for w are all ones
%  x is 2D
%  y is 1D
%  w is 1D
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

global sizeImage;
% it contains four elements [1d,2d,c]

c=sizeImage(3);
% the number of coils

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

if(~isfield(opts,'w'))
    w=ones(n,1);
else
    w=opts.w;
    if norm( length(w) - n ) ~=0
        error('\n opts.w is incorrect');
    end
end


if(isfield(opts,'x0'))
    x0=opts.x0;
else
    x0=zeros(n,c);
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
    
    if (opts.pFlag==1)
        temp=sqrt(sum(abs(ATy).^2,2) );
        
        lambda_max=max(temp./w);
    end
    
    if (opts.pFlag==2)
        for i=1:size(ATy,2)
            V(:,i)=opts.W( ATy(:,i) );
        end
        temp=sqrt(sum(abs(V).^2,2) );
        
        lambda_max=max(temp);
    end
            
    if (opts.pFlag==3)
        for i=1:size(ATy,2)
            V(:,i)=opts.W( ATy(:,i) );
        end
        temp=abs(V);
        
        lambda_max=max(temp,[],1);
    end
else
    if ( (opts.pFlag==1) || (opts.pFlag==2) )
        lambda_max=1;
    elseif (opts.pFlag==3)
        lambda_max=ones(1,c);
    end
end

%% VV
for i=1:size(ATy,2)
    V(:,i)=opts.W( ATy(:,i) );
end

if(~isfield(opts,'pFlag'))
    opts.pFlag=1; % 1 for regular (L2,1) norm 
                  % 2 for redundant Harr
                  % 3 for independent 
end

if (opts.pFlag==1)
    proximity_operator = @(x, w, th) (soft_2D(x, w, th));
elseif (opts.pFlag==2)
    % proximity_operator = @(x, w, th) operator_rh(x, 1, th, opts.W, opts.WT);
    iter=50;
    p=zeros(size(V));
    q=p;
    
    %proximity_operator = @(v, w, th) Dykstra_sparse_group(opts.W, opts.WT, v,th,iter);
    proximity_operator = @(v, w, th) Dykstra_group(opts.W, opts.WT, v,th,iter);
elseif (opts.pFlag==3)
    proximity_operator = @(x, w, th) operator_rh2(x, 1, th, opts.W, opts.WT);
end


fprintf('\n lambda_max=%e',lambda_max);

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
                
        lambda0 = 0.1*lambda_max;
        
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
% a 1D vector

% assign xp with x, and Axp with Ax
xp=x; Axp=Ax; xxp=zeros(n,c);

% alphap and alpha are used for computing the weight in forming search point
alphap=0; alpha=1;

for iterStep=1:opts.maxIter
    fprintf('\n fista_slep_2D_sos iteration number: %d', iterStep);
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
        
        x = proximity_operator(v, w, lambda/L);
        
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
    
%     figure(1);
%     imshow( reshape( sqrt(sum(abs(x).^2,2)) ,[256,256]) , []);
%     pause(0.01);
    
    % --------------------------- step 3 ---------------------------
    % update alpha and alphap, and check whether converge
    alphap=alpha; alpha= (1+ sqrt(4*alpha*alpha +1))/2;
    
    xxp=x-xp;   Axy=Ax-y;
    
    if(nargout>=2)
        if ( (opts.pFlag==1) || (opts.pFlag==2) )
            funVal(iterStep)=norm(Axy(:))^2/2 + ...
                lambda * sum( w.* sqrt( sum(abs(x).^2,2) ) );       
        end
        if ( opts.pFlag==3 )
            funVal(iterStep)=norm(Axy(:))^2/2;   
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
                    funVal(iterStep)=norm(Axy(:))^2/2 + ...
                        lambda * sum( w.* sqrt( sum(abs(x).^2,2) ) );
                    %fprintf('iterStep: %d, lambda =%e, function value = %e', iterStep, lambda, funVal(iterStep));
                end
                
                if (abs( funVal(iterStep) - funVal(iterStep-1) ) <= opts.tol)
                    break;
                end
            end
        case 1
            if iterStep>=2
                
                if(nargout<2)
                    funVal(iterStep)=norm(Axy(:))^2/2 + ...
                        lambda * sum( w.* sqrt( sum(abs(x).^2,2) ) );
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
    function x=Dykstra_group(W, WT, v,th,iter)
        % solving the problem
        %    0.5 \|x-v\|_2^2 + lambda \|W x\|_1
        % via Dykstra's algorithm
        

        %fprintf('\n lambda=%e',lambda);
        
        for i=1:size(v,2)
            Wv(:,i)=W( v(:,i) );
        end
        
        xx=Wv-p-q;        
        
        if iterStep<20
            iter=1;
        else
            iter=10;
        end
        
        for i=1:iter
            
            yy=soft_2D(xx+p,1, th);
            
            p=p+(xx-yy);
            
            
            yyq=yy+q;
            
            for j=1:size(yyq,2)
                xx(:,j)=W( WT( yyq(:,j) ));
            end
                        
            q=q+(yy-xx);
            
            %fprintf('\n i=%d, %e, %e',i, norm(xx(:)-yy(:)), norm(xx(:)));
            
            if norm(xx-yy)< 1e-4
                break;
            end
        end
        
        for i=1:size(xx,2)
            x(:,i)=WT( xx(:,i) );
        end
    end


    function x=Dykstra_sparse_group(W, WT, v,th,iter)
        % solving the problem
        %    0.5 \|x-v\|_2^2 + lambda \|W x\|_1
        % via Dykstra's algorithm
        
        
        fprintf('\n lambda=%e',lambda);
        
        for i=1:size(v,2)
            Wv(:,i)=W( v(:,i) );
        end
        
        xx=Wv-p-q;
        
        if iterStep<20
            iter=1;
        else
            iter=10;
        end
        
        for i=1:iter
            
            yy=soft_2D( soft(xx+p,th) ,1, th);
            
            p=p+(xx-yy);
            
            
            yyq=yy+q;
            
            for j=1:size(yyq,2)
                xx(:,j)=W( WT( yyq(:,j) ));
            end
            
            q=q+(yy-xx);
            
            %fprintf('\n i=%d, %e, %e',i, norm(xx(:)-yy(:)), norm(xx(:)));
            
            if norm(xx-yy)< 1e-4
                break;
            end
        end
        
        for i=1:size(xx,2)
            x(:,i)=WT( xx(:,i) );
        end
    end
end
