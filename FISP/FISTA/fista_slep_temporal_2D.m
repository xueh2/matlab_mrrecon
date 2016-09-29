function [x, funVal,  ValueL, mse]=fista_slep_temporal_2D(A, AT, y, W, WT, lambda, opts)
%
%%
% Function LeastR
%      Least Squares with Nesterov's Approach
%
%% Problem
%
%  min  1/2 \|Ax -y\|^2 + lambda \|Wx\|_1 
%
%% Input parameters:
%
%  A, AT, W, WT-         operator
%  y,c -      vector (of size nx1)
%
%% Output parameters:
%  x-         Solution
%  funVal-    Function value during iterations
%
%%
% Composed by Jun Liu on March 22, 2011

%% Verify and initialize the parameters
%%

% Verify the number of input parameters
if (nargin <7)
    error('\n Inputs: A, AT, y, W, WT, lambda, opts should be specified!\n');
end

if(isfield(opts,'n'))
    n=opts.n;
else
    error('opts.n should be specified');
end


global sizeReconImage;
global sizeImage;

if(~isfield(opts,'pFlag'))
    opts.pFlag=1;
end

%opts.pFlag
% 1: redundant Harr
% 2: TV
% 3: pure L1
% 4: redundant Haar (treat real and imag separately)


% ATy
ATy=AT(y);

% WATy
WATy=W(ATy);

if(isfield(opts,'rFlag'))
    if (opts.rFlag~=0 && opts.rFlag~=1)
        error('opts.rFlag error');
    end
else
    opts.rFlag=0; % no lambda_max
end

if (opts.rFlag)
    
    if (opts.pFlag==1 || opts.pFlag==4)
        lambda_max=max(abs(WATy));
        lambda_max2=max(abs(ATy));
        
        %[lambda_max2,lambda_max/2, lambda_max2- lambda_max/2, 0.99*lambda_max2- lambda_max/2]
    else
        lambda_max=max(abs(ATy));
        lambda_max2=lambda_max;
    end
else
    lambda_max=1;
end


scale=opts.scale;


lambda_max

lambda=lambda*lambda_max;

if(isfield(opts,'x0'))
    x0=opts.x0;
else
    if (isfield(opts,'init'))
        if (opts.init==0)
            x0=zeros(n,1);
        else
            % we set x0=tau * ATy 
            tau=max(0, (2* norm(ATy)^2 - lambda * sum(abs(WATy)))/ 2 / norm(A(ATy))^2 );
            x0=tau*ATy;
        end
    else        
        x0=zeros(n,1);
    end
end

if(isfield(opts,'L'))
    L=opts.L;
else
    L=0.1; % We assume that the maximum eigenvalue of A'A is over 1
end

if(~isfield(opts,'maxIter'))
    opts.maxIter=1000;
end

if(~isfield(opts,'tFlag'))
    opts.tFlag=3;
end

if(~isfield(opts,'tol'))
    opts.tol=1e-12;
end

if (opts.pFlag==4)
    if(isfield(opts,'d'))
        d=opts.d; 
    else
        d=ones(size(WATy));
    end
end


if ((opts.pFlag==1) || (opts.pFlag==4) )
    %proximity_operator = @(x, th) WT(soft(W(x), th));
    p=zeros(length(WATy),1);
    q=p;
    proximity_operator = @(x,th) WT( Dykstra(W(x),th));
end

if (opts.pFlag==3)
    proximity_operator = @(x, th) (soft(x, th));
end


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
        
        if(isfield(opts,'maxFlag'))
            if opts.maxFlag==1
                norminf=lambda_max;
            else
                norminf=lambda_max2;
            end
        else
            norminf=lambda_max2;
        end
        
        lambda0 = 0.9*norminf;
        
        if(isfield(opts,'lambda_bar'))
            lambda_bar=opts.lambda_bar *norminf;
        else            
            lambda_bar = 0.001*norminf;
        end
        
        L=0.1;
        
        eta=0.8;
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
Ax=A(x);

% assign xp with x, and Axp with Ax
xp=x; Axp=Ax; xxp=zeros(n,1);

% alphap and alpha are used for computing the weight in forming search point
alphap=0; alpha=1;

mm=sizeReconImage(1);
nn=sizeReconImage(2);
Z0r=zeros(mm ,nn);
Z0c=zeros(mm,nn);

for iterStep=1:opts.maxIter
    %fprintf('\n fista_slep iteration number: %d', iterStep);
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
        
        if (opts.pFlag==1 || opts.pFlag==3)
            x = proximity_operator(v, lambda/L);
        elseif(opts.pFlag==2)
            V=reshape(abs(v), sizeReconImage);
            
            Wr=[ones(mm-1,nn)*lambda/L; zeros(1,nn)];
            Wc=[ones(mm,nn-1)*lambda/L, zeros(mm,1)]';
            
            [X,Zr,Zc]=mat_tv2(V,Wr,Wc,Z0r,Z0c,30);
            x=X(:);
            
            Z0r=Zr; Z0c=Zc;
        end
        
        if (opts.pFlag==4)
            x = proximity_operator(v, lambda/L*d);
        end
        
        v=x-s;  % the difference between the new approximate solution x
        % and the search point s
        
        % compute A x
        Ax=A(x);
        
        Av=Ax -As;
        r_sum=norm(v)^2; l_sum=norm(Av)^2;
        
        if (r_sum <=1e-20) && (iterStep>=20)
            bFlag=1; % this shows that, the gradient step makes little improvement
            break;
        end
        
        % the condition is ||Av||_2^2 <= (L) * ||v||_2^2
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
        if (opts.pFlag==1)
            funVal(iterStep)=norm(Axy)^2/2 +lambda * sum(abs(W(x)));
        elseif (opts.pFlag==3)
            funVal(iterStep)=norm(Axy)^2/2 +lambda * sum(abs(x));
        else
            funVal(iterStep)=norm(Axy)^2/2;
        end
        
        %fprintf('iterStep: %d, lambda =%e, function value = %e', iterStep, lambda, funVal(iterStep));
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
    
    if (iterStep>=20)
        
        switch(opts.tFlag)
            case 0
                if iterStep>=2
                    if (abs( funVal(iterStep) - funVal(iterStep-1) ) <= opts.tol)
                        break;
                    end
                end
            case 1
                if iterStep>=2
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
                norm_xxp=norm(xxp);
                if ( norm_xxp <=opts.tol)
                    break;
                end
            case 4
                norm_xp=norm(xp);    norm_xxp=norm(xxp);
                if ( norm_xxp <=opts.tol * max(norm_xp,1))
                    break;
                end
            case 5
                if iterStep>=opts.maxIter
                    break;
                end
        end
    end
end

function x=Dykstra(Wv,th)

%     p=zeros(length(Wv),1);
%     q=p;
%     x=Wv;
    x=Wv-p-q;
    

    if (lambda > lambda_bar )
        p=0;
        q=0;
        x=Wv;
    else
        p=sign(p).*min(abs(p),th);
        x=Wv-p-q;
    end
    
    for i=1:10
        xplusp=reshape( x+p, 2*[sizeImage(1), sizeImage(2)]);
        yy=soft_temporal_2d(xplusp,th, scale);
        
        p=p+(x-yy);
        
        x=W( WT(yy+q) );
        q=q+(yy-x);
        
        %fprintf('\n %e', 0.5* norm(x-Wv).^2 + th * sum(abs(x)));
        
        if norm(x-yy)< norm(x)* 1e-4
            break;
        end
        
    end
end

end
