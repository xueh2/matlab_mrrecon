function mjelp
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ============================================================
%       Purpose: This program computes Jacobian elliptic functions
%                sn u, cn u and dn u using subroutine JELP
%       Input  : u   --- Argument of Jacobian elliptic fuctions
%                Hk  --- Modulus k(0 ó k ó 1)
%       Output : ESN --- sn u
%                ECN --- cn u
%                EDN --- dn u
%                EPH --- phi(in degrees)
%       Example:
%                k = .5,(K(k)= 1.68575035), and u = u0*K
%                u0       phi       sn u        cn u        dn u
%              ----------------------------------------------------
%               0.0      .0000    .0000000   1.0000000   1.0000000
%               0.5    47.0586    .7320508    .6812500    .9306049
%               1.0    90.0000   1.0000000    .0000000    .8660254
%               1.5   132.9414    .7320508   -.6812500    .9306049
%               2.0   180.0000    .0000000  -1.0000000   1.0000000
%       ============================================================
u=[];hk=[];esn=[];ecn=[];edn=[];eph=[];
fprintf(1,'%s \n','please enter k and u ');
%        READ(*,*)HK,U
hk=.5;
u=1.0.*1.68575035;
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   k        u          phi        sn u');fprintf(1,'%s \n', '        cn u        dn u');
fprintf(1,'%s ',' -------------------------------------');fprintf(1,'%s \n','---------------------------');
[u,hk,esn,ecn,edn,eph]=jelp(u,hk,esn,ecn,edn,eph);
fprintf(1,[repmat(' ',1,1),'%5.3g','%12.7g',repmat(' ',1,2),'%9.5g',repmat('%12.7g',1,3) ' \n'],hk,u,eph,esn,ecn,edn);
%format(1x,f5.3,f12.7,2x,f9.5,3f12.7);
end
function [u,hk,esn,ecn,edn,eph]=jelp(u,hk,esn,ecn,edn,eph,varargin);
%       ========================================================
%       Purpose: Compute Jacobian elliptic functions sn u, cn u
%                and dn u
%       Input  : u   --- Argument of Jacobian elliptic fuctions
%                Hk  --- Modulus k(0 ó k ó 1)
%       Output : ESN --- sn u
%                ECN --- cn u
%                EDN --- dn u
%                EPH --- phi(in degrees)
%       ========================================================
 r=zeros(1,40);
pi=3.14159265358979d0;
a0=1.0d0;
b0=sqrt(1.0d0-hk.*hk);
for  n=1:40;
a=(a0+b0)./2.0d0;
b=sqrt(a0.*b0);
c=(a0-b0)./2.0d0;
r(n)=c./a;
if(c < 1.0d-7)break; end;
a0=a;
b0=b;
end;
dn=2.0d0.^n.*a.*u;
for  j=n:-1:1;
t=r(j).*sin(dn);
sa=atan(t./sqrt(abs(1.0d0-t.*t)));
d=.5d0.*(dn+sa);
dn=d;
end;  j=1-1;
eph=d.*180.0d0./pi;
esn=sin(d);
ecn=cos(d);
edn=sqrt(1.0d0-hk.*hk.*esn.*esn);
return;
end

