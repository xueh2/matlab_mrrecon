function mcerzo
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     ===============================================================
%     Purpose : This program evaluates the complex zeros of error
%     function erf(z)using subroutine CERZO
%     Input:    NT --- Total number of zeros
%     Example:  NT = 10
%     n     complex zeros of erf(z)n     complex zeros of erf(z)
%     ------------------------------------------------------------------
%     -
%     1   1.450616163 + i 1.880943000   6   4.158998400 + i 4.435571444
%     2   2.244659274 + i 2.616575141   7   4.516319400 + i 4.780447644
%     3   2.839741047 + i 3.175628100   8   4.847970309 + i 5.101588043
%     4   3.335460735 + i 3.646174376   9   5.158767908 + i 5.403332643
%     5   3.769005567 + i 4.060697234  10   5.452192201 + i 5.688837437
%     ===============================================================
nt=[];zo=[];
 zo=zeros(1,100);
fprintf(1,'%s \n','please enter nt ');
%      READ(*,*)NT
nt=10;
fprintf(1,[repmat(' ',1,2),'nt=','%3g' ' \n'],nt);
[nt,zo]=cerzo(nt,zo);
fprintf(1,'%s \n','  *****    please wait ...    *****');
fprintf(1,'%0.15g \n');
fprintf(1,'%s \n','  n        complex zeros of erf(z)');
fprintf(1,'%s \n','-------------------------------------');
for  i=1:nt;
fprintf(1,[repmat(' ',1,1),'%3g',repmat(' ',1,2),'%13.8g',repmat(' ',1,2),'+i','%13.8g' ' \n'],i,zo(i));
end;  i=nt+1;
%format(2x,',i3);
%format(1x,i3,2x,f13.8,2x,2h+i,f13.8);
end
function [nt,zo]=cerzo(nt,zo,varargin);
%     ===============================================================
%     Purpose : Evaluate the complex zeros of error function erf(z)
%     using the modified Newton's iteration method
%     Input :   NT --- Total number of zeros
%     Output:   ZO(L)--- L-th zero of erf(z), L=1,2,...,NT
%     Routine called: CERF for computing erf(z)and erf'(z)
%     ===============================================================
z=[];zf=[];zd=[];
zo(:)=0.0;
w=0.0;
pi=3.141592653589793d0;
for  nr=1:nt;
pu=sqrt(pi.*(4.0d0.*nr-0.5d0));
pv=pi.*sqrt(2.0d0.*nr-0.25d0);
px=0.5.*pu-0.5.*log(pv)./pu;
py=0.5.*pu+0.5.*log(pv)./pu;
z=complex(px,py);
it=0;
while(1==1);
it=it+1;
[z,zf,zd]=cerf(z,zf,zd);
zp=complex(1.0d0,0.0d0);
for  i=1:nr-1;
zp=zp.*(z-zo(i));
end;  i=nr-1+1;
zfd=zf./zp;
zq=complex(0.0d0,0.0d0);
for  i=1:nr-1;
zw=complex(1.0d0,0.0d0);
for  j=1:nr-1;
if(j ~= i);
zw=zw.*(z-zo(j));
end;
end;  j=nr-1+1;
zq=zq+zw;
end;  i=nr-1+1;
zgd=(zd-zq.*zfd)./zp;
z=z-zfd./zgd;
w0=w;
w=abs(z);
if(it <= 50&abs((w-w0)./w)> 1.0d-11)break; end;
end;
zo(nr)=z;
end;
return;
end
function [z,cer,cder]=cerf(z,cer,cder,varargin);
%     ==========================================================
%     Purpose: Compute complex Error function erf(z)& erf'(z)
%     Input:   z   --- Complex argument of erf(z)
%     x   --- Real part of z
%     y   --- Imaginary part of z
%     Output:  CER --- erf(z)
%     CDER --- erf'(z)
%     ==========================================================
w=0.0;
w1=0.0;
w2=0.0;
eps=1.0d-12;
pi=3.141592653589793d0;
x=real(z);
y=imag(z);
x2=x.*x;
if(x <= 3.5d0);
er=1.0d0;
r=1.0d0;
for  k=1:100;
r=r.*x2./(k+0.5d0);
er=er+r;
if(abs(er-w)<= eps.*abs(er))break; end;
w=er;
end;
c0=2.0d0./sqrt(pi).*x.*exp(-x2);
er0=c0.*er;
else;
er=1.0d0;
r=1.0d0;
for  k=1:12;
r=-r.*(k-0.5d0)./x2;
er=er+r;
end;  k=12+1;
c0=exp(-x2)./(x.*sqrt(pi));
er0=1.0d0-c0.*er;
end;
if(y == 0.0d0);
err=er0;
eri=0.0d0;
else;
cs=cos(2.0d0.*x.*y);
ss=sin(2.0d0.*x.*y);
er1=exp(-x2).*(1.0d0-cs)./(2.0d0.*pi.*x);
ei1=exp(-x2).*ss./(2.0d0.*pi.*x);
er2=0.0d0;
for  n=1:100;
er2=er2+exp(-.25d0.*n.*n)./(n.*n+4.0d0.*x2).*(2.0d0.*x-2.0d0.*x.*cosh(n.*y).*cs+n.*sinh(n.*y).*ss);
if(abs((er2-w1)./er2)< eps)break; end;
w1=er2;
end;
c0=2.0d0.*exp(-x2)./pi;
err=er0+er1+c0.*er2;
ei2=0.0d0;
for  n=1:100;
ei2=ei2+exp(-.25d0.*n.*n)./(n.*n+4.0d0.*x2).*(2.0d0.*x.*cosh(n.*y).*ss+n.*sinh(n.*y).*cs);
if(abs((ei2-w2)./ei2)< eps)break; end;
w2=ei2;
end;
eri=ei1+c0.*ei2;
end;
cer=complex(err,eri);
cder=2.0d0./sqrt(pi).*exp(-z.*z);
return;
end

