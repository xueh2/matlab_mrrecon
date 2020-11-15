function mairyb
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ============================================================
%       Purpose: This program computes Airy functions and their
%                derivatives using subroutine AIRYB
%       Input:   x  --- Argument of Airy function
%       Output:  AI --- Ai(x)
%                BI --- Bi(x)
%                AD --- Ai'(x)
%                BD --- Bi'(x)
%       Example:
%   x       Ai(x)Bi(x)Ai'(x)Bi'(x)
%  ----------------------------------------------------------------
%   0   .35502805D+00  .61492663D+00 -.25881940D+00  .44828836D+00
%  10   .11047533D-09  .45564115D+09 -.35206337D-09  .14292361D+10
%  20   .16916729D-26  .21037650D+26 -.75863916D-26  .93818393D+26
%  30   .32082176D-48  .90572885D+47 -.17598766D-47  .49533045D+48
%   x       Ai(-x)Bi(-x)Ai'(-x)Bi'(-x)
%  ----------------------------------------------------------------
%   0       .35502805      .61492663     -.25881940      .44828836
%  10       .04024124     -.31467983      .99626504      .11941411
%  20      -.17640613     -.20013931      .89286286     -.79142903
%  30      -.08796819     -.22444694     1.22862060     -.48369473
%       ============================================================
x=[];ai=[];bi=[];ad=[];bd=[];
fprintf(1,'%s \n','please enter x ');
%        READ(*,*)X
x=20;
[x,ai,bi,ad,bd]=airyb(x,ai,bi,ad,bd);
fprintf(1,[repmat(' ',1,4),'x',repmat(' ',1,8),'ai(x)',repmat(' ',1,11),'bi(x)',repmat(' ',1,11),'ai''(x)',repmat(' ',1,10),'bi''(x)' ' \n']);
fprintf(1,[repmat(' ',1,2),'----------------------------------','-----------------------------------' ' \n']);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%16.8g',1,4) ' \n'],x,ai,bi,ad,bd);
fprintf(1,'%0.15g \n');
[dumvar1,ai,bi,ad,bd]=airyb(-x,ai,bi,ad,bd);
fprintf(1,[repmat(' ',1,4),'x',repmat(' ',1,8),'ai(-x)',repmat(' ',1,10),'bi(-x)',repmat(' ',1,10),'ai''(-x)',repmat(' ',1,9),'bi''(-x)' ' \n']);
fprintf(1,[repmat(' ',1,2),'----------------------------------','-----------------------------------' ' \n']);
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%16.8g',1,4) ' \n'],x,ai,bi,ad,bd);
%format(1x,f5.1,4d16.8);
%format(1x,f5.1,4d16.8);
%format(4x,'x',8x,'ai(x)',11x,'bi(x)',11x,'ai''(x)', 10x,'bi''(x)');
%format(2x,'----------------------------------','-----------------------------------');
%format(4x,'x',8x,'ai(-x)',10x,'bi(-x)',10x, 'ai''(-x)',9x,'bi''(-x)');
end
function [x,ai,bi,ad,bd]=airyb(x,ai,bi,ad,bd,varargin);
%       =======================================================
%       Purpose: Compute Airy functions and their derivatives
%       Input:   x  --- Argument of Airy function
%       Output:  AI --- Ai(x)
%                BI --- Bi(x)
%                AD --- Ai'(x)
%                BD --- Bi'(x)
%       =======================================================
 ck=zeros(1,41);
dk=zeros(1,41);
eps=1.0d-15;
pi=3.141592653589793d0;
c1=0.355028053887817d0;
c2=0.258819403792807d0;
sr3=1.732050807568877d0;
xa=abs(x);
xq=sqrt(xa);
if(x > 0.0d0)xm=5.0; end;
if(x <= 0.0d0)xm=8.0; end;
if(x == 0.0d0);
ai=c1;
bi=sr3.*c1;
ad=-c2;
bd=sr3.*c2;
return;
end;
if(xa <= xm);
fx=1.0d0;
r=1.0d0;
for  k=1:40;
r=r.*x./(3.0d0.*k).*x./(3.0d0.*k-1.0d0).*x;
fx=fx+r;
if(abs(r)< abs(fx).*eps)break; end;
end;
gx=x;
r=x;
for  k=1:40;
r=r.*x./(3.0d0.*k).*x./(3.0d0.*k+1.0d0).*x;
gx=gx+r;
if(abs(r)< abs(gx).*eps)break; end;
end;
ai=c1.*fx-c2.*gx;
bi=sr3.*(c1.*fx+c2.*gx);
df=0.5d0.*x.*x;
r=df;
for  k=1:40;
r=r.*x./(3.0d0.*k).*x./(3.0d0.*k+2.0d0).*x;
df=df+r;
if(abs(r)< abs(df).*eps)break; end;
end;
dg=1.0d0;
r=1.0d0;
for  k=1:40;
r=r.*x./(3.0d0.*k).*x./(3.0d0.*k-2.0d0).*x;
dg=dg+r;
if(abs(r)< abs(dg).*eps)break; end;
end;
ad=c1.*df-c2.*dg;
bd=sr3.*(c1.*df+c2.*dg);
else;
xe=xa.*xq./1.5d0;
xr1=1.0d0./xe;
xar=1.0d0./xq;
xf=sqrt(xar);
rp=0.5641895835477563d0;
r=1.0d0;
for  k=1:40;
r=r.*(6.0d0.*k-1.0d0)./216.0d0.*(6.0d0.*k-3.0d0)./k.*(6.0d0.*k-5.0d0)./(2.0d0.*k-1.0d0);
ck(k)=r;
dk(k)=-(6.0d0.*k+1.0d0)./(6.0d0.*k-1.0d0).*ck(k);
end;  k=40+1;
km=fix(24.5-xa);
if(xa < 6.0)km=14; end;
if(xa > 15.0)km=10; end;
if(x > 0.0d0);
sai=1.0d0;
sad=1.0d0;
r=1.0d0;
for  k=1:km;
r=-r.*xr1;
sai=sai+ck(k).*r;
sad=sad+dk(k).*r;
end;  k=km+1;
sbi=1.0d0;
sbd=1.0d0;
r=1.0d0;
for  k=1:km;
r=r.*xr1;
sbi=sbi+ck(k).*r;
sbd=sbd+dk(k).*r;
end;  k=km+1;
xp1=exp(-xe);
ai=0.5d0.*rp.*xf.*xp1.*sai;
bi=rp.*xf./xp1.*sbi;
ad=-.5d0.*rp./xf.*xp1.*sad;
bd=rp./xf./xp1.*sbd;
else;
xcs=cos(xe+pi./4.0d0);
xss=sin(xe+pi./4.0d0);
ssa=1.0d0;
sda=1.0d0;
r=1.0d0;
xr2=1.0d0./(xe.*xe);
for  k=1:km;
r=-r.*xr2;
ssa=ssa+ck(2.*k).*r;
sda=sda+dk(2.*k).*r;
end;  k=km+1;
ssb=ck(1).*xr1;
sdb=dk(1).*xr1;
r=xr1;
for  k=1:km;
r=-r.*xr2;
ssb=ssb+ck(2.*k+1).*r;
sdb=sdb+dk(2.*k+1).*r;
end;  k=km+1;
ai=rp.*xf.*(xss.*ssa-xcs.*ssb);
bi=rp.*xf.*(xcs.*ssa+xss.*ssb);
ad=-rp./xf.*(xcs.*sda+xss.*sdb);
bd=rp./xf.*(xss.*sda-xcs.*sdb);
end;
end;
return;
end

