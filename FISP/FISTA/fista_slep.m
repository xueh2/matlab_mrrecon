function [x, funVal,  ValueL, mse]=fista_slep(A, AT, y, W, WT, lambda, opts)
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

if(~isfield(opts,'pFlag'))
    opts.pFlag=1;
end

%opts.pFlag
% 1: redundant Harr
% 2: TV
% 3: pure L1


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
    
    if (opts.pFlag==1)
        lambda_max=max(abs(WATy));
        lambda_max2=max(abs(ATy));
        
        [lambda_max2,lambda_max]
    else
        lambda_max=max(abs(ATy));
        lambda_max2=lambda_max;
    end
else
    lambda_max=1;
end

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
    L=1; % We assume that the maximum eigenvalue of A'A is over 1
end

if(~isfield(opts,'maxIter'))
    opts.maxIter=1000;
end

if(~isfield(opts,'tFlag'))
    opts.tFlag=3;
end

if(~isfield(opts,'tol'))
    opts.tol=1e-6;
end



if (opts.pFlag==1)
    proximity_operator = @(x, th) WT(soft(W(x), th));
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
        
        lambda0 = 1.5*norminf;
        
        if(isfield(opts,'lambda_bar'))
            lambda_bar=opts.lambda_bar *norminf;
        else            
            lambda_bar = 0.001*norminf;
        end
        
        L=0.1;
        
        eta=0.8;
        lambda = lambda0;
    end
else
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
    fprintf('\n fista_slep iteration number: %d', iterStep);
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
        
        v=x-s;  % the difference between the new approximate solution x
        % and the search point s
        
        % compute A x
        Ax=A(x);
        
        Av=Ax -As;
        r_sum=norm(v)^2; l_sum=norm(Av)^2;
        
        if (r_sum <=1e-20)
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





% %% adaptive line search
% 
% % .mFlag=1, and .lFlag=1
% %  refomulate the problem as the constrained convex optimization
% %  problem, and then apply adaptive line search scheme
% 
% % Problem:
% %    min  1/2 || A x - y||^2 + 1/2 rsL2 * ||x||_2^2 + z * t' * 1
% %    s.t.   |x| <= t
% 
% if(opts.mFlag==1 && opts.lFlag==1)
%     
%     bFlag=0; % this flag tests whether the gradient step only changes a little
% 
%     L=1 + rsL2;
%     % We assume that the maximum eigenvalue of A'A is over 1
%     
%     gamma=1;
%     % we shall set the value of gamma = L,
%     % where L is appropriate for the starting point
% 
%     xp=x; Axp=Ax;
%     % store x and Ax
%     xxp=zeros(n,1);
%     % the difference of x and xp    
%     t=abs(x); tp=t;
%     % t is the upper bound of absolute value of x
%     
%     % compute AT Ax
%     if (opts.nFlag==0)
%         ATAx=A'*Ax;
%     elseif (opts.nFlag==1)
%         ATAx=A'*Ax - sum(Ax) * mu';  ATAx=ATAx./nu;
%     else
%         invNu=Ax./nu;                ATAx=A'*invNu-sum(invNu)*mu';
%     end
%     
%     % We begin the adaptive line search in the following
%     %
%     % Note that, in the line search, L and beta are changing
%     
%     for iterStep=1:opts.maxIter
% 
%         ATAxp=ATAx;
%         % store ATAx to ATAxp
% 
%         if (iterStep~=1)
%             % compute AT Ax
%             if (opts.nFlag==0)
%                 ATAx=A'*Ax;
%             elseif (opts.nFlag==1)
%                 ATAx=A'*Ax - sum(Ax) * mu';  ATAx=ATAx./nu;
%             else
%                 invNu=Ax./nu;                ATAx=A'*invNu-sum(invNu)*mu';
%             end
%         end
% 
%         %--------- Line Search for L begins
%         while (1)
%             if (iterStep~=1)
%                 alpha= (-gamma+ sqrt(gamma*gamma + 4* L * gamma)) / (2*L);
%                 beta= (gamma - gamma* alphap) / (alphap * gamma + alphap* L * alpha);
%                 % beta is the coefficient for generating search point s
% 
%                 s=x + beta* xxp;   s_t= t + beta * (t -tp);
%                 As=Ax + beta* (Ax-Axp);
%                 ATAs=ATAx + beta * (ATAx- ATAxp);
%                 % compute the search point s, A * s, and A' * A * s
%             else
%                 alpha= (-1+ sqrt(5)) / 2;
%                 beta=0; s=x; s_t=t; As=Ax; ATAs=ATAx;
%             end
% 
%             g=ATAs-ATy + rsL2 * s;
%             % compute the gradient g
%            
%             % let s walk in a step in the antigradient of s 
%             u=s-g/L;
%             v= s_t - lambda / L;
% 
%             % projection
%             [xnew, tnew]=ep1R(u,v,n);
% 
%             v=xnew-s;  % the difference between the new approximate solution x
%                             % and the search point s
%             v_t=tnew-s_t;
%             
%             % compute A xnew
%             if (opts.nFlag==0)
%                 Axnew=A* xnew;
%             elseif (opts.nFlag==1)
%                 invNu=xnew./nu; mu_invNu=mu * invNu;
%                 Axnew=A*invNu -repmat(mu_invNu, m, 1);
%             else
%                 Axnew=A*xnew-repmat(mu*xnew, m, 1);     Axnew=Axnew./nu;
%             end
% 
%             Av=Axnew -As;
%             r_sum=v'*v + v_t'*v_t; l_sum=Av'*Av + v'*v * rsL2;
%             
%             if (r_sum <=1e-20)
%                 bFlag=1; % this shows that, the gradient step makes little improvement
%                 break;
%             end
%             
%             % the condition is ||Av||_2^2 + rsL2 * ||v||_2^2
%             %                       <= L * (||v||_2^2 + ||v_t|| _2^2 )
%             if(l_sum <= r_sum * L)
%                 break;
%             else
%                 L=max(2*L, l_sum/r_sum);
%                 % fprintf('\n L=%5.6f',L);
%             end
%         end
%         %--------- Line Search for L ends
% 
%         gamma=L* alpha* alpha;    alphap=alpha;
%         % update gamma, and alphap
%         
%         ValueL(iterStep)=L;
% 
%         tao=L * r_sum / l_sum;
%         if (tao >=5)
%             L=L*0.8;
%         end
%         % decrease the value of L
% 
%         xp=x;    x=xnew; xxp=x-xp;
%         Axp=Ax;  Ax=Axnew;
%         % update x and Ax with xnew and Axnew        
%         tp=t; t=tnew;
%         % update tp and t       
%         
%         Axy=Ax-y;
%         funVal(iterStep)=Axy' * Axy/2 + rsL2/2 * x'*x + lambda * sum(t);
%         % compute function value
%         
%         if (bFlag)
%             % fprintf('\n The program terminates as the gradient step changes the solution very small.');
%             break;
%         end
% 
%         switch(opts.tFlag)
%             case 0
%                 if iterStep>=2
%                     if (abs( funVal(iterStep) - funVal(iterStep-1) ) <= opts.tol)
%                         break;
%                     end
%                 end
%             case 1
%                 if iterStep>=2
%                     if (abs( funVal(iterStep) - funVal(iterStep-1) ) <=...
%                             opts.tol* funVal(iterStep-1))
%                         break;
%                     end
%                 end
%             case 2
%                 if ( funVal(iterStep)<= opts.tol)
%                     break;
%                 end
%             case 3
%                 norm_xxp=sqrt(xxp'*xxp+ norm(t-tp)^2);
%                 if ( norm_xxp <=opts.tol)
%                     break;
%                 end
%             case 4
%                 norm_xp=sqrt(xp'*xp + tp'*tp);    norm_xxp=sqrt(xxp'*xxp+ norm(t-tp)^2);
%                 if ( norm_xxp <=opts.tol * max(norm_xp,1))
%                     break;
%                 end
%             case 5
%                 if iterStep>=opts.maxIter
%                     break;
%                 end
%         end
%     end
% end
% 
% 
% %%
% if(opts.mFlag==0 && opts.lFlag==1)
%     error('\n The function does not support opts.mFlag=0 & opts.lFlag=1!');
% end