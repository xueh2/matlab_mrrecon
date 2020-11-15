function mairyzo
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =========================================================
%       Purpose: This program computes the first NT zeros of Airy
%                functions Ai(x)and Ai'(x), and the associated
%                values of Ai(a')and Ai'(a), and the first NT
%                zeros of Airy functions Bi(x)and Bi'(x), and
%                the associated values of Bi(b')and Bi'(b)using
%                subroutine AIRYZO
%       Input :  NT    --- Total number of zeros
%                KF    --- Function code
%                          KF=1 for Ai(x)and Ai'(x)
%                          KF=2 for Bi(x)and Bi'(x)
%       Output:  XA(m)--- a, the m-th zero of Ai(x)or
%                          b, the m-th zero of Bi(x)
%                XB(m)--- a', the m-th zero of Ai'(x)or
%                          b', the m-th zero of Bi'(x)
%                XC(m)--- Ai(a')or Bi(b')
%                XD(m)--- Ai'(a)or Bi'(b)
%(m --- Serial number of zeros)
%       Example: NT=5
%       m         a            Ai'(a)a'          Ai(a')
%      -----------------------------------------------------------
%       1    -2.33810741     .70121082   -1.01879297    .53565666
%       2    -4.08794944    -.80311137   -3.24819758   -.41901548
%       3    -5.52055983     .86520403   -4.82009921    .38040647
%       4    -6.78670809    -.91085074   -6.16330736   -.35790794
%       5    -7.94413359     .94733571   -7.37217726    .34230124
%       m         b            Bi'(b)b'          Bi(b')
%      -----------------------------------------------------------
%       1    -1.17371322     .60195789   -2.29443968   -.45494438
%       2    -3.27109330    -.76031014   -4.07315509    .39652284
%       3    -4.83073784     .83699101   -5.51239573   -.36796916
%       4    -6.16985213    -.88947990   -6.78129445    .34949912
%       5    -7.37676208     .92998364   -7.94017869   -.33602624
%       ==========================================================
nt=[];kf=[];xa=[];xb=[];xc=[];xd=[];
 xa=zeros(1,50);
xb=zeros(1,50);
xc=zeros(1,50);
xd=zeros(1,50);
fprintf(1,[repmat(' ',1,10),'kf=1 for ai(x)and ai''(x); kf=2 for bi(x)',' and bi''(x)' ' \n']);
fprintf(1,[repmat(' ',1,10),'nt is the number of the zeros' ' \n']);
fprintf(1,'%s \n','please enter kf,nt ');
%        READ(*,*)KF,NT
kf=1;
nt=5;
fprintf(1,[repmat(' ',1,1),'kf=','%2g',',     ','nt=','%3g' ' \n'],kf,nt);
if(kf == 1);
fprintf(1,'%s ','  m        a             ai''(a)a''');fprintf(1,'%s \n','           ai(a'')');
elseif(kf == 2);
fprintf(1,'%s ','  m        b             bi''(b)b''');fprintf(1,'%s \n','           bi(b'')');
end;
fprintf(1,'%s ','---------------------------------');fprintf(1,'%s \n','---------------------------');
[nt,kf,xa,xb,xc,xd]=airyzo(nt,kf,xa,xb,xc,xd);
for  k=1:nt;
fprintf(1,[repmat(' ',1,1),'%3g',repmat(' ',1,1),repmat('%14.8g',1,3),'%13.8g' ' \n'],k,xa(k),xd(k),xb(k),xc(k));
end;  k=nt+1;
%format(1x,i3,1x,3f14.8,f13.8);
%format(1x,,i2,',     ',,i3);
%format(10x,'kf=1 for ai(x)and ai''(x); kf=2 for bi(x)',' and bi''(x)');
%format(10x,'nt is the number of the zeros');
end
function [nt,kf,xa,xb,xc,xd]=airyzo(nt,kf,xa,xb,xc,xd,varargin);
%       ========================================================
%       Purpose: Compute the first NT zeros of Airy functions
%                Ai(x)and Ai'(x), a and a', and the associated
%                values of Ai(a')and Ai'(a); and the first NT
%                zeros of Airy functions Bi(x)and Bi'(x), b and
%                b', and the associated values of Bi(b')and
%                Bi'(b)
%       Input :  NT    --- Total number of zeros
%                KF    --- Function code
%                          KF=1 for Ai(x)and Ai'(x)
%                          KF=2 for Bi(x)and Bi'(x)
%       Output:  XA(m)--- a, the m-th zero of Ai(x)or
%                          b, the m-th zero of Bi(x)
%                XB(m)--- a', the m-th zero of Ai'(x)or
%                          b', the m-th zero of Bi'(x)
%                XC(m)--- Ai(a')or Bi(b')
%                XD(m)--- Ai'(a)or Bi'(b)
%(m --- Serial number of zeros)
%       Routine called: AIRYB for computing Airy functions and
%                       their derivatives
%       =======================================================
x=[];ai=[];bi=[];ad=[];bd=[];
pi=3.141592653589793d0;
for  i=1:nt;
if(kf == 1);
u=3.0.*pi.*(4.0.*i-1)./8.0d0;
u1=1./(u.*u);
rt0=-(u.*u).^(1.0./3.0).*((((-15.5902.*u1+.929844).*u1-.138889).*u1+.10416667d0).*u1+1.0d0);
elseif(kf == 2);
if(i == 1);
rt0=-1.17371;
else;
u=3.0.*pi.*(4.0.*i-3.0)./8.0;
u1=1.0d0./(u.*u);
rt0=-(u.*u).^(1.0./3.0).*((((-15.5902.*u1+.929844).*u1-.138889).*u1+.10416667).*u1+1.0);
end;
end;
rt=1.0e300;
while(abs((rt-rt0)./rt)> 1.d-9);
x=rt0;
[x,ai,bi,ad,bd]=airyb(x,ai,bi,ad,bd);
if(kf == 1)rt=rt0-ai./ad; end;
if(kf == 2)rt=rt0-bi./bd; end;
if(abs((rt-rt0)./rt)> 1.d-9);
rt0=rt;
end;
end;
xa(i)=rt;
if(kf == 1)xd(i)=ad; end;
if(kf == 2)xd(i)=bd; end;
end;  i=fix(nt)+1;
for  i=1:nt;
if(kf == 1);
if(i == 1);
rt0=-1.01879;
else;
u=3.0.*pi.*(4.0.*i-3.0)./8.0;
u1=1./(u.*u);
rt0=-(u.*u).^(1.0./3.0).*((((15.0168.*u1-.873954).*u1+.121528).*u1-.145833d0).*u1+1.0d0);
end;
elseif(kf == 2);
if(i == 1);
rt0=-2.29444;
else;
u=3.0.*pi.*(4.0.*i-1.0)./8.0;
u1=1.0./(u.*u);
rt0=-(u.*u).^(1.0./3.0).*((((15.0168.*u1-.873954).*u1+.121528).*u1-.145833).*u1+1.0);
end;
end;
rt=1.0e300;
while(abs((rt-rt0)./rt)> 1.0d-9);
x=rt0;
[x,ai,bi,ad,bd]=airyb(x,ai,bi,ad,bd);
if(kf == 1)rt=rt0-ad./(ai.*x); end;
if(kf == 2)rt=rt0-bd./(bi.*x); end;
if(abs((rt-rt0)./rt)> 1.0d-9);
rt0=rt;
end;
end;
xb(i)=rt;
if(kf == 1)xc(i)=ai; end;
if(kf == 2)xc(i)=bi; end;
end;  i=fix(nt)+1;
return;
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
if(abs(r./fx)< eps)break; end;
end;
gx=x;
r=x;
for  k=1:40;
r=r.*x./(3.0d0.*k).*x./(3.0d0.*k+1.0d0).*x;
gx=gx+r;
if(abs(r./gx)< eps)break; end;
end;
ai=c1.*fx-c2.*gx;
bi=sr3.*(c1.*fx+c2.*gx);
df=.5d0.*x.*x;
r=df;
for  k=1:40;
r=r.*x./(3.0d0.*k).*x./(3.0d0.*k+2.0d0).*x;
df=df+r;
if(abs(r./df)< eps)break; end;
end;
dg=1.0d0;
r=1.0d0;
for  k=1:40;
r=r.*x./(3.0d0.*k).*x./(3.0d0.*k-2.0d0).*x;
dg=dg+r;
if(abs(r./dg)< eps)break; end;
end;
ad=c1.*df-c2.*dg;
bd=sr3.*(c1.*df+c2.*dg);
else;
xe=xa.*xq./1.5d0;
xr1=1.0d0./xe;
xar=1.0d0./xq;
xf=sqrt(xar);
rp=.5641895835477563d0;
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
ai=.5d0.*rp.*xf.*xp1.*sai;
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

