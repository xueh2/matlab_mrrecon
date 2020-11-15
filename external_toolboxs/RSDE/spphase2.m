%  spphase2 
%
% The code is obtained at: http://dollar.biz.uiowa.edu/col/ye/matlab.html
%
%  This is the Phase 2 procedure called by SPSOLQP. 

 lamda=(1.-beta)*lamda;
 go=0;
 dx = ones(n,1)./x;
%
%  Repeatly solve an ellipsoid constrained QP problem by solving a linear
%  system equation until find a positive solution.
%
 while go <= 0,
   DD = sparse(1:n,1:n,(lamda*dx).*dx,n,n);
%
   u=[Q+DD A';A sparse(m,m)]\[-(Q*x+c);sparse(m,1)];
   xx=x+u(1:n);
   go=min(xx);
   if go > 0,
     ob=xx'*Q*xx/2+c'*xx;
     go = min([go obvalue-ob+eps]);
   end;
   lamda=2*lamda;
   if lamda >= (1+abs(obvalue))/toler,
     disp('The problem seems unbounded.');
     return
   end; 
 end;
%
 y=-u(n+1:n+m);
 u=u(1:n);
 nora = min(u./x);
 if nora < 0,
   nora=-alpha/nora;
 elseif nora == 0,
   nora=alpha;
 else
   nora=inf;
 end
%
 w1 = u'*Q*u;
 w2 = -u'*(Q*x+c);
 if w1 > 0,
  nora=min([w2/w1,nora]);
 end;
 if nora == inf,
  ob = -inf;
 else
   x =x+nora*u;
   ob=x'*Q*x/2+c'*x;
 end;
 clear u dx xx DD w1 w2



