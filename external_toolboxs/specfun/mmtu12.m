function mmtu12
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%
%
%     ===============================================================
%     Purpose: This program computes the modified Mathieu functions
%     of the first and second kinds, Mcm(1,2,x,q)and
%     Msm(1,2,x,q), and their derivatives using
%     subroutine MTU12
%     Input:   KF --- Function code
%     KF=1 for computing Mcm(x,q)
%     KF=2 for computing Msm(x,q)
%     KC --- Function Code
%     KC=1 for computing Mcm(1,x,q)and Mcm(1)'(x,q)
%     or Msm(1,x,q)and Msm(1)'(x,q)
%     KC=2 for computing Mcm(2,x,q)and Mcm(2)'(x,q)
%     or Msm(2,x,q)and Msm(2)'(x,q)
%     KC=3 for both modified Mathieu functions of the
%     first and second kinds, and their
%     derivatives
%     m  --- Order of Mathieu functions
%     q  --- Parameter of Mathieu functions
%     x  --- Argument of Mathieu functions
%     Output:  F1R --- Mcm(1,x,q)or Msm(1,x,q)
%     D1R --- Derivative of Mcm(1,x,q)or Msm(1,x,q)
%     F2R --- Mcm(2,x,q)or Msm(2,x,q)
%     D2R --- Derivative of Mcm(2,x,q)or Msm(2,x,q)
%     ===============================================================
%
%
kf=[];kc=[];m=[];q=[];x=[];f1r=[];d1r=[];f2r=[];d2r=[];
 d1r=0;
 d2r=0;
 f1r=0;
 f2r=0;
 q=0;
 x=0;
 kc=0;
 kf=0;
 m=0;
fprintf(1,'%s \n','please enter kf, m, q and x ');
%     READ(*,*)KF,M,Q,X
kf=1;
m=1;
q=2;
x=3;
fprintf(1,[repmat(' ',1,1),'kf =','%2g',',  ','m =','%3g',',  ','q =','%5.1g',',  ','x =','%5.1g' ' \n'],fix(kf) , fix(m) , q , x);
%format(1x,',i2,',  ',',i3,',  ',',f5.1,',  ',',f5.1);
kc=3;
[kf,kc,m,q,x,f1r,d1r,f2r,d2r]=mtu12(fix(kf),fix(kc),fix(m),q,x,f1r,d1r,f2r,d2r);
fprintf(1,'%0.15g \n');
if(kf==1);
fprintf(1,'%s ','   x      mcm(1,x,q)mcm(1)''(x,q)' );fprintf(1,'%s \n','    mcm(2,x,q)mcm(2)''(x,q)');
else;
fprintf(1,'%s ','   x      msm(1,x,q)msm(1)''(x,q)' );fprintf(1,'%s \n','    msm(2,x,q)msm(2)''(x,q)');
end;
fprintf(1,'%s ',' --------------------------------------' );fprintf(1,'%s \n','-------------------------------');
fprintf(1,[repmat(' ',1,1),'%5.1g',repmat('%16.8g',1,4) ' \n'],x , f1r , d1r , f2r , d2r);
%format(1x,f5.1,4d16.8);
fprintf(1,'%0.15g \n');
fprintf(1,[repmat(' ',1,1),'wronskian=','%16.8g',repmat(' ',1,3),'should equal   2/pi=','%16.8g' ' \n'],f1r.*d2r-f2r.*d1r , .63661977236758d0);
%format(1x,',e16.8,3x,',e16.8);
fprintf(1,[repmat(' ',1,1), '\n ' ,repmat(' ',1,1),'caution: this check is not accurate if it ','involves', '\n ' ,repmat(' ',1,1),'         the subtraction of two ','similar numbers' ' \n']);
%format(1x,./1x,['caution: this check is not accurate if it '],'involves',./1x,'         the subtraction of two ','similar numbers');
end
%--/mtu12.f90  processed by SPAG 6.53Rc at 13:11 on 28 Jun 2004
%
%
function [kf,kc,m,q,x,f1r,d1r,f2r,d2r]=mtu12(kf,kc,m,q,x,f1r,d1r,f2r,d2r,varargin);
%
%     ==============================================================
%     Purpose: Compute modified Mathieu functions of the first and
%     second kinds, Mcm(1,2,x,q)and Msm(1,2,x,q),
%     and their derivatives
%     Input:   KF --- Function code
%     KF=1 for computing Mcm(x,q)
%     KF=2 for computing Msm(x,q)
%     KC --- Function Code
%     KC=1 for computing the first kind
%     KC=2 for computing the second kind
%     or Msm(2,x,q)and Msm(2)'(x,q)
%     KC=3 for computing both the first
%     and second kinds
%     m  --- Order of Mathieu functions
%     q  --- Parameter of Mathieu functions(q ò 0)
%     x  --- Argument of Mathieu functions
%     Output:  F1R --- Mcm(1,x,q)or Msm(1,x,q)
%     D1R --- Derivative of Mcm(1,x,q)or Msm(1,x,q)
%     F2R --- Mcm(2,x,q)or Msm(2,x,q)
%     D2R --- Derivative of Mcm(2,x,q)or Msm(2,x,q)
%     Routines called:
%(1)CVA2 for computing the characteristic values
%(2)FCOEF for computing expansion coefficients
%(3)JYNB for computing Jn(x), Yn(x)and their
%     derivatives
%     ==============================================================
%
%
kd=[];a=[];fg=[];km=[];u1=[];nm=[];bj1=[];dj1=[];by1=[];dy1=[];u2=[];bj2=[];dj2=[];by2=[];dy2=[];
 a=0;
 bj1=zeros(1,251+1);
 bj2=zeros(1,251+1);
 by1=zeros(1,251+1);
 by2=zeros(1,251+1);
c1=0;
 c2=0;
 dj1=zeros(1,251+1);
 dj2=zeros(1,251+1);
 dy1=zeros(1,251+1);
dy2=zeros(1,251+1);
 eps=0;
 fg=zeros(1,251);
 qm=0;
 u1=0;
 u2=0;
 w1=0;
 w2=0;
 ic=0;
 k=0;
 kd=0;
 km=0;
 nm=0;
eps=1.0d-14;
if(kf==1 & m==2.*fix(m./2));
kd=1;
end;
if(kf==1 & m~=2.*fix(m./2));
kd=2;
end;
if(kf==2 & m~=2.*fix(m./2));
kd=3;
end;
if(kf==2 & m==2.*fix(m./2));
kd=4;
end;
[kd,m,q,a]=cva2(fix(kd),fix(m),q,a);
if(q<=1.0d0);
qm=7.5+56.1.*sqrt(q)-134.7.*q+90.7.*sqrt(q).*q;
else;
qm=17.0+3.1.*sqrt(q)-.126.*q+.0037.*sqrt(q).*q;
end;
km=fix(fix(qm+0.5.*fix(m)));
[kd,m,q,a,fg]=fcoef(fix(kd),fix(m),q,a,fg);
ic=fix(fix(fix(m)./2)+1);
if(kd==4);
ic=fix(fix(m)./2);
end;
c1=exp(-x);
c2=exp(x);
u1=sqrt(q).*c1;
u2=sqrt(q).*c2;
[km,u1,nm,bj1,dj1,by1,dy1]=jynb(fix(km),u1,fix(nm),bj1,dj1,by1,dy1);
[km,u2,nm,bj2,dj2,by2,dy2]=jynb(fix(km),u2,fix(nm),bj2,dj2,by2,dy2);
if(kc~=2);
f1r=0.0d0;
for k=1 : km;
if(kd==1);
f1r=f1r+(-1).^(fix(ic)+fix(k)).*fg(fix(k)).*bj1(fix(k)-1+1).*bj2(fix(k)-1+1);
elseif(kd==2 | kd==3);
f1r=f1r+(-1).^(fix(ic)+fix(k)).*fg(fix(k)).*(bj1(fix(k)-1+1).*bj2(fix(k)+1)+(-1).^fix(kd).*bj1(fix(k)+1).*bj2(fix(k)-1+1));
else;
f1r=f1r+(-1).^(fix(ic)+fix(k)).*fg(fix(k)).*(bj1(fix(k)-1+1).*bj2(fix(k)+1+1)-bj1(fix(k)+1+1).*bj2(fix(k)-1+1));
end;
if(k>=5 & abs(f1r-w1)<abs(f1r).*eps);
break;
end;
w1=f1r;
end;
f1r=f1r./fg(1);
d1r=0.0d0;
for k=1 : km;
if(kd==1);
d1r=d1r+(-1).^(fix(ic)+fix(k)).*fg(fix(k)).*(c2.*bj1(fix(k)-1+1).*dj2(fix(k)-1+1)-c1.*dj1(fix(k)-1+1).*bj2(fix(k)-1+1));
elseif(kd==2 | kd==3);
d1r=d1r+(-1).^(fix(ic)+fix(k)).*fg(fix(k)).*(c2.*(bj1(fix(k)-1+1).*dj2(fix(k)+1)+(-1).^fix(kd).*bj1(fix(k)+1).*dj2(fix(k)-1+1))-c1.*(dj1(fix(k)-1+1).*bj2(fix(k)+1)+(-1).^fix(kd).*dj1(fix(k)+1).*bj2(fix(k)-1+1)));
else;
d1r=d1r+(-1).^(fix(ic)+fix(k)).*fg(fix(k)).*(c2.*(bj1(fix(k)-1+1).*dj2(fix(k)+1+1)-bj1(fix(k)+1+1).*dj2(fix(k)-1+1))-c1.*(dj1(fix(k)-1+1).*bj2(fix(k)+1+1)-dj1(fix(k)+1+1).*bj2(fix(k)-1+1)));
end;
if(k>=5 & abs(d1r-w2)<abs(d1r).*eps);
break;
end;
w2=d1r;
end;
d1r=d1r.*sqrt(q)./fg(1);
if(kc==1);
return;
end;
end;
f2r=0.0d0;
for k=1 : km;
if(kd==1);
f2r=f2r+(-1).^(fix(ic)+fix(k)).*fg(fix(k)).*bj1(fix(k)-1+1).*by2(fix(k)-1+1);
elseif(kd==2 | kd==3);
f2r=f2r+(-1).^(fix(ic)+fix(k)).*fg(fix(k)).*(bj1(fix(k)-1+1).*by2(fix(k)+1)+(-1).^fix(kd).*bj1(fix(k)+1).*by2(fix(k)-1+1));
else;
f2r=f2r+(-1).^(fix(ic)+fix(k)).*fg(fix(k)).*(bj1(fix(k)-1+1).*by2(fix(k)+1+1)-bj1(fix(k)+1+1).*by2(fix(k)-1+1));
end;
if(k>=5 & abs(f2r-w1)<abs(f2r).*eps);
break;
end;
w1=f2r;
end;
f2r=f2r./fg(1);
d2r=0.0d0;
for k=1 : km;
if(kd==1);
d2r=d2r+(-1).^(fix(ic)+fix(k)).*fg(fix(k)).*(c2.*bj1(fix(k)-1+1).*dy2(fix(k)-1+1)-c1.*dj1(fix(k)-1+1).*by2(fix(k)-1+1));
elseif(kd==2 | kd==3);
d2r=d2r+(-1).^(fix(ic)+fix(k)).*fg(fix(k)).*(c2.*(bj1(fix(k)-1+1).*dy2(fix(k)+1)+(-1).^fix(kd).*bj1(fix(k)+1).*dy2(fix(k)-1+1))-c1.*(dj1(fix(k)-1+1).*by2(fix(k)+1)+(-1).^fix(kd).*dj1(fix(k)+1).*by2(fix(k)-1+1)));
else;
d2r=d2r+(-1).^(fix(ic)+fix(k)).*fg(fix(k)).*(c2.*(bj1(fix(k)-1+1).*dy2(fix(k)+1+1)-bj1(fix(k)+1+1).*dy2(fix(k)-1+1))-c1.*(dj1(fix(k)-1+1).*by2(fix(k)+1+1)-dj1(fix(k)+1+1).*by2(fix(k)-1+1)));
end;
if(k>=5 & abs(d2r-w2)<abs(d2r).*eps);
break;
end;
w2=d2r;
end;
d2r=d2r.*sqrt(q)./fg(1);
end
%--/fcoef.f90  processed by SPAG 6.53Rc at 13:11 on 28 Jun 2004
%
%
function [kd,m,q,a,fc]=fcoef(kd,m,q,a,fc,varargin);
%
%     =====================================================
%     Purpose: Compute expansion coefficients for Mathieu
%     functions and modified Mathieu functions
%     Input :  m  --- Order of Mathieu functions
%     q  --- Parameter of Mathieu functions
%     KD --- Case code
%     KD=1 for cem(x,q,m = 0,2,4,...)
%     KD=2 for cem(x,q,m = 1,3,5,...)
%     KD=3 for sem(x,q,m = 1,3,5,...)
%     KD=4 for sem(x,q,m = 2,4,6,...)
%     A  --- Characteristic value of Mathieu
%     functions for given m and q
%     Output:  FC(k)--- Expansion coefficients of Mathieu
%     functions(k= 1,2,...,KM)
%     FC(1),FC(2),FC(3),... correspond to
%     A0,A2,A4,... for KD=1 case, A1,A3,
%     A5,... for KD=2 case, B1,B3,B5,...
%     for KD=3 case and B2,B4,B6,... for
%     KD=4 case
%     =====================================================
%
%
 f=0;
 f1=0;
 f2=0;
 f3=0;
 qm=0;
 s=0;
 s0=0;
 sp=0;
 ss=0;
 u=0;
 v=0;
 i=0;
 ifoo=0;
 ifoo2=0;
 j=0;
 k=0;
 kb=0;
 km=0;
for i=1 : 251;
fc(i)=0.0d0;
end; i=fix( 251+1);
if(q<=1.0d0);
qm=7.5+56.1.*sqrt(q)-134.7.*q+90.7.*sqrt(q).*q;
else;
qm=17.0+3.1.*sqrt(q)-.126.*q+.0037.*sqrt(q).*q;
end;
km=fix(fix(qm+0.5.*fix(m)));
for k=1 : km+1;
fc(k)=0.0d0;
end; k=fix( fix(km)+1+1);
if(q==0.0d0);
if(kd==1);
fc((m+2)./2)=1.0d0;
if(m==0);
fc(1)=1.0d0./sqrt(2.0d0);
end;
elseif(kd==4);
fc(m./2)=1.0d0;
else;
fc((m+1)./2)=1.0d0;
end;
return;
end;
kb=0;
s=0.0d0;
f=1.0d-100;
u=0.0d0;
fc(km)=0.0d0;
ifoo=1;
while(ifoo==1);
if(kd==1);
for k=km : -1: 3 ;
v=u;
u=f;
f=(a-4.0d0.*fix(k).*fix(k)).*u./q-v;
if(abs(f)<abs(fc(k+1)));
kb=fix(fix(k));
fc(1)=1.0d-100;
sp=0.0d0;
f3=fc(fix(k)+1);
fc(2)=a./q.*fc(1);
fc(3)=(a-4.0d0).*fc(2)./q-2.0d0.*fc(1);
u=fc(2);
f1=fc(3);
for i=3 : kb;
v=u;
u=f1;
f1=(a-4.0d0.*(fix(i)-1.0d0).^2).*u./q-v;
fc(i+1)=f1;
if(i==kb);
f2=f1;
end;
if(i~=kb);
sp=sp+f1.*f1;
end;
end; i=fix( fix(kb)+1);
sp=sp+2.0d0.*fc(1).^2+fc(2).^2+fc(3).^2;
ss=s+sp.*(f3./f2).^2;
s0=sqrt(1.0d0./ss);
for j=1 : km;
if(j<=kb+1);
fc(j)=s0.*fc(fix(j)).*f3./f2;
else;
fc(j)=s0.*fc(fix(j));
end;
end; j=fix( fix(km)+1);
ifoo=0;
break;
%GO TO 85
else;
fc(k)=f;
s=s+f.*f;
end;
end;
if(ifoo==0);
break;
end;
fc(2)=q.*fc(3)./(a-4.0d0-2.0d0.*q.*q./a);
fc(1)=q./a.*fc(2);
s=s+2.0d0.*fc(1).^2+fc(2).^2;
s0=sqrt(1.0d0./s);
for k=1 : km;
fc(k)=s0.*fc(fix(k));
end; k=fix( fix(km)+1);
elseif(kd==2 | kd==3);
ifoo2=1;
for k=km : -1: 3 ;
v=u;
u=f;
f=(a-(2.0d0.*fix(k)-1).^2).*u./q-v;
if(abs(f)>=abs(fc(k)));
fc(k-1)=f;
s=s+f.*f;
else;
kb=fix(fix(k));
f3=fc(fix(k));
ifoo2=0;
break;
end;
end;
if(ifoo2==1);
fc(1)=q./(a-1.0d0-(-1).^fix(kd).*q).*fc(2);
s=s+fc(1).*fc(1);
s0=sqrt(1.0d0./s);
for k=1 : km;
fc(k)=s0.*fc(fix(k));
end; k=fix( fix(km)+1);
break;
%GO TO 85
end;
fc(1)=1.0d-100;
fc(2)=(a-1.0d0-(-1).^fix(kd).*q)./q.*fc(1);
sp=0.0d0;
u=fc(1);
f1=fc(2);
for i=2 : kb-1;
v=u;
u=f1;
f1=(a-(2.0d0.*fix(i)-1.0d0).^2).*u./q-v;
if(i~=kb-1);
fc(i+1)=f1;
sp=sp+f1.*f1;
else;
f2=f1;
end;
end; i=fix( fix(kb)-1+1);
sp=sp+fc(1).^2+fc(2).^2;
ss=s+sp.*(f3./f2).^2;
s0=1.0d0./sqrt(ss);
for j=1 : km;
if(j<kb);
fc(j)=s0.*fc(fix(j)).*f3./f2;
end;
if(j>=kb);
fc(j)=s0.*fc(fix(j));
end;
end; j=fix( fix(km)+1);
elseif(kd==4);
ifoo2=1;
for k=km : -1: 3 ;
v=u;
u=f;
f=(a-4.0d0.*fix(k).*fix(k)).*u./q-v;
if(abs(f)>=abs(fc(k)));
fc(k-1)=f;
s=s+f.*f;
else;
kb=fix(fix(k));
f3=fc(fix(k));
ifoo2=0;
break;
end;
end;
if(ifoo2==1);
fc(1)=q./(a-4.0d0).*fc(2);
s=s+fc(1).*fc(1);
s0=sqrt(1.0d0./s);
for k=1 : km;
fc(k)=s0.*fc(fix(k));
end; k=fix( fix(km)+1);
break;
%GO TO 85
end;
fc(1)=1.0d-100;
fc(2)=(a-4.0d0)./q.*fc(1);
sp=0.0d0;
u=fc(1);
f1=fc(2);
for i=2 : kb-1;
v=u;
u=f1;
f1=(a-4.0d0.*fix(i).*fix(i)).*u./q-v;
if(i~=kb-1);
fc(i+1)=f1;
sp=sp+f1.*f1;
else;
f2=f1;
end;
end; i=fix( fix(kb)-1+1);
sp=sp+fc(1).^2+fc(2).^2;
ss=s+sp.*(f3./f2).^2;
s0=1.0d0./sqrt(ss);
for j=1 : km;
if(j<kb);
fc(j)=s0.*fc(fix(j)).*f3./f2;
end;
if(j>=kb);
fc(j)=s0.*fc(fix(j));
end;
end; j=fix( fix(km)+1);
end;
ifoo=0;
end;
if(fc(1)<0.0d0);
for j=1 : km;
fc(j)=-fc(fix(j));
end; j=fix( fix(km)+1);
end;
end
%--/cva2.f90  processed by SPAG 6.53Rc at 13:11 on 28 Jun 2004
%
%
function [kd,m,q,a]=cva2(kd,m,q,a,varargin);
%
%     ======================================================
%     Purpose: Calculate a specific characteristic value of
%     Mathieu functions
%     Input :  m  --- Order of Mathieu functions
%     q  --- Parameter of Mathieu functions
%     KD --- Case code
%     KD=1 for cem(x,q,m = 0,2,4,...)
%     KD=2 for cem(x,q,m = 1,3,5,...)
%     KD=3 for sem(x,q,m = 1,3,5,...)
%     KD=4 for sem(x,q,m = 2,4,6,...)
%     Output:  A  --- Characteristic value
%     Routines called:
%(1)REFINE for finding accurate characteristic
%     value using an iteration method
%(2)CV0 for finding initial characteristic
%     values using polynomial approximation
%(3)CVQM for computing initial characteristic
%     values for q ó 3*m
%(3)CVQL for computing initial characteristic
%     values for q ò m*m
%     ======================================================
%
%
q1=[];a1=[];q2=[];a2=[];qq=[];iflag=[];
 a1=0;
 a2=0;
 delta=0;
 q1=0;
 q2=0;
 qq=0;
 i=0;
 iflag=0;
 ndiv=0;
 nn=0;
if(m<=12 | q<=3.0.*m | q>m.*m);
[kd,m,q,a]=cv0(fix(kd),fix(m),q,a);
if(q~=0.0d0);
[kd,m,q,a]=refine(fix(kd),fix(m),q,a,1);
end;
else;
ndiv=10;
delta=(fix(m)-3.0).*fix(m)./fix(ndiv);
if((q-3.0.*m)<=(m.*m-q));
while (1);
nn=fix(fix((q-3.0.*fix(m))./delta)+1);
delta=(q-3.0.*fix(m))./fix(nn);
q1=2.0.*fix(m);
[m,q1,a1]=cvqm(fix(m),q1,a1);
q2=3.0.*fix(m);
[m,q2,a2]=cvqm(fix(m),q2,a2);
qq=3.0.*fix(m);
for i=1 : nn;
qq=qq+delta;
a=(a1.*q2-a2.*q1+(a2-a1).*qq)./(q2-q1);
iflag=1;
if(i==nn);
iflag=-1;
end;
[kd,m,qq,a,iflag]=refine(fix(kd),fix(m),qq,a,fix(iflag));
q1=q2;
q2=qq;
a1=a2;
a2=a;
end; i=fix( fix(nn)+1);
if(iflag==-10);
ndiv=fix(fix(ndiv).*2);
delta=(fix(m)-3.0).*fix(m)./fix(ndiv);
break;
end;
end;
else;
while (1);
nn=fix(fix((fix(m).*fix(m)-q)./delta)+1);
delta=(fix(m).*fix(m)-q)./fix(nn);
q1=fix(m).*(fix(m)-1.0);
[kd,m,q1,a1]=cvql(fix(kd),fix(m),q1,a1);
q2=fix(m).*fix(m);
[kd,m,q2,a2]=cvql(fix(kd),fix(m),q2,a2);
qq=fix(m).*fix(m);
for i=1 : nn;
qq=qq-delta;
a=(a1.*q2-a2.*q1+(a2-a1).*qq)./(q2-q1);
iflag=1;
if(i==nn);
iflag=-1;
end;
[kd,m,qq,a,iflag]=refine(fix(kd),fix(m),qq,a,fix(iflag));
q1=q2;
q2=qq;
a1=a2;
a2=a;
end; i=fix( fix(nn)+1);
if(iflag==-10);
ndiv=fix(fix(ndiv).*2);
delta=(fix(m)-3.0).*fix(m)./fix(ndiv);
break;
end;
end;
end;
end;
end
%--/refine.f90  processed by SPAG 6.53Rc at 13:11 on 28 Jun 2004
%
%
function [kd,m,q,a,iflag]=refine(kd,m,q,a,iflag,varargin);
%
%     =====================================================
%     Purpose: calculate the accurate characteristic value
%     by the secant method
%     Input :  m --- Order of Mathieu functions
%     q --- Parameter of Mathieu functions
%     A --- Initial characteristic value
%     Output:  A --- Refineed characteristic value
%     Routine called:  CVF for computing the value of F for
%     characteristic equation
%     ========================================================
%
%
x0=[];mj=[];f0=[];x1=[];f1=[];x=[];f=[];
 eps=0;
 f=0;
 f0=0;
 f1=0;
 x=0;
 x0=0;
 x1=0;
 it=0;
 mj=0;
eps=1.0d-14;
mj=fix(10+fix(m));
x0=a;
[kd,m,q,x0,mj,f0]=cvf(fix(kd),fix(m),q,x0,fix(mj),f0);
x1=1.002.*a;
[kd,m,q,x1,mj,f1]=cvf(fix(kd),fix(m),q,x1,fix(mj),f1);
for it=1 : 100;
mj=fix(fix(mj)+1);
x=x1-(x1-x0)./(1.0d0-f0./f1);
[kd,m,q,x,mj,f]=cvf(fix(kd),fix(m),q,x,fix(mj),f);
if(abs(1.0-x1./x)<eps | f==0.0);
break;
end;
x0=x1;
f0=f1;
x1=x;
f1=f;
end;
a=x;
end
%--/cvf.f90  processed by SPAG 6.53Rc at 13:11 on 28 Jun 2004
%
%
function [kd,m,q,a,mj,f]=cvf(kd,m,q,a,mj,f,varargin);
%
%     ======================================================
%     Purpose: Compute the value of F for characteristic
%     equation of Mathieu functions
%     Input :  m --- Order of Mathieu functions
%     q --- Parameter of Mathieu functions
%     A --- Characteristic value
%     Output:  F --- Value of F for characteristic equation
%     ======================================================
%
%
 b=0;
 t0=0;
 t1=0;
 t2=0;
 ic=0;
 j=0;
 j0=0;
 jf=0;
 l=0;
 l0=0;
b=a;
ic=fix(fix(fix(m)./2));
l=0;
l0=0;
j0=2;
jf=fix(fix(ic));
if(kd==1);
l0=2;
end;
if(kd==1);
j0=3;
end;
if(kd==2 | kd==3);
l=1;
end;
if(kd==4);
jf=fix(fix(ic)-1);
end;
t1=0.0d0;
for j=mj : -1: ic+1 ;
t1=-q.*q./((2.0d0.*fix(j)+fix(l)).^2-b+t1);
end; j=fix( fix(ic)+1 -1);
if(m<=2);
t2=0.0d0;
if(kd==1 & m==0);
t1=t1+t1;
end;
if(kd==1 & m==2);
t1=-2.0.*q.*q./(4.0-b+t1)-4.0;
end;
if(kd==2 & m==1);
t1=t1+q;
end;
if(kd==3 & m==1);
t1=t1-q;
end;
else;
if(kd==1);
t0=4.0d0-b+2.0d0.*q.*q./b;
end;
if(kd==2);
t0=1.0d0-b+q;
end;
if(kd==3);
t0=1.0d0-b-q;
end;
if(kd==4);
t0=4.0d0-b;
end;
t2=-q.*q./t0;
for j=j0 : jf;
t2=-q.*q./((2.0d0.*fix(j)-fix(l)-fix(l0)).^2-b+t2);
end; j=fix( fix(jf)+1);
end;
f=(2.0d0.*fix(ic)+fix(l)).^2+t1+t2-b;
end
%--/cv0.f90  processed by SPAG 6.53Rc at 13:11 on 28 Jun 2004
%
%
function [kd,m,q,a0]=cv0(kd,m,q,a0,varargin);
%
%     =====================================================
%     Purpose: Compute the initial characteristic value of
%     Mathieu functions for m ó 12  or q ó 300 or
%     q ò m*m
%     Input :  m  --- Order of Mathieu functions
%     q  --- Parameter of Mathieu functions
%     Output:  A0 --- Characteristic value
%     Routines called:
%(1)CVQM for computing initial characteristic
%     value for q ó 3*m
%(2)CVQL for computing initial characteristic
%     value for q ò m*m
%     ====================================================
%
%
 q2=0;
q2=q.*q;
if(m==0);
if(q<=1.0);
a0=(((.0036392.*q2-.0125868).*q2+.0546875).*q2-.5).*q2;
elseif(q<=10.0);
a0=((3.999267d-3.*q-9.638957d-2).*q-.88297).*q+.5542818;
else;
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
end;
elseif(m==1);
if(q<=1.0 & kd==2);
a0=(((-6.51e-4.*q-.015625).*q-.125).*q+1.0).*q+1.0;
elseif(q<=1.0 & kd==3);
a0=(((-6.51e-4.*q+.015625).*q-.125).*q-1.0).*q+1.0;
elseif(q<=10.0 & kd==2);
a0=(((-4.94603d-4.*q+1.92917d-2).*q-.3089229).*q+1.33372).*q+.811752;
elseif(q<=10.0 & kd==3);
a0=((1.971096d-3.*q-5.482465d-2).*q-1.152218).*q+1.10427;
else;
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
end;
elseif(m==2);
if(q<=1.0 & kd==1);
a0=(((-.0036391.*q2+.0125888).*q2-.0551939).*q2+.416667).*q2+4.0;
elseif(q<=1.0 & kd==4);
a0=(.0003617.*q2-.0833333).*q2+4.0;
elseif(q<=15 & kd==1);
a0=(((3.200972d-4.*q-8.667445d-3).*q-1.829032d-4).*q+.9919999).*q+3.3290504;
elseif(q<=10.0 & kd==4);
a0=((2.38446d-3.*q-.08725329).*q-4.732542d-3).*q+4.00909;
else;
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
end;
elseif(m==3);
if(q<=1.0 & kd==2);
a0=((6.348e-4.*q+.015625).*q+.0625).*q2+9.0;
elseif(q<=1.0 & kd==3);
a0=((6.348e-4.*q-.015625).*q+.0625).*q2+9.0;
elseif(q<=20.0 & kd==2);
a0=(((3.035731d-4.*q-1.453021d-2).*q+.19069602).*q-.1039356).*q+8.9449274;
elseif(q<=15.0 & kd==3);
a0=((9.369364d-5.*q-.03569325).*q+.2689874).*q+8.771735;
else;
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
end;
elseif(m==4);
if(q<=1.0 & kd==1);
a0=((-2.1e-6.*q2+5.012e-4).*q2+.0333333).*q2+16.0;
elseif(q<=1.0 & kd==4);
a0=((3.7e-6.*q2-3.669e-4).*q2+.0333333).*q2+16.0;
elseif(q<=25.0 & kd==1);
a0=(((1.076676d-4.*q-7.9684875d-3).*q+.17344854).*q-.5924058).*q+16.620847;
elseif(q<=20.0 & kd==4);
a0=((-7.08719d-4.*q+3.8216144d-3).*q+.1907493).*q+15.744;
else;
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
end;
elseif(m==5);
if(q<=1.0 & kd==2);
a0=((6.8e-6.*q+1.42e-5).*q2+.0208333).*q2+25.0;
elseif(q<=1.0 & kd==3);
a0=((-6.8e-6.*q+1.42e-5).*q2+.0208333).*q2+25.0;
elseif(q<=35.0 & kd==2);
a0=(((2.238231d-5.*q-2.983416d-3).*q+.10706975).*q-.600205).*q+25.93515;
elseif(q<=25.0 & kd==3);
a0=((-7.425364d-4.*q+2.18225d-2).*q+4.16399d-2).*q+24.897;
else;
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
end;
elseif(m==6);
if(q<=1.0);
a0=(.4d-6.*q2+.0142857).*q2+36.0;
elseif(q<=40.0 & kd==1);
a0=(((-1.66846d-5.*q+4.80263d-4).*q+2.53998d-2).*q-.181233).*q+36.423;
elseif(q<=35.0 & kd==4);
a0=((-4.57146d-4.*q+2.16609d-2).*q-2.349616d-2).*q+35.99251;
else;
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
end;
elseif(m==7);
if(q<=10.0);
[m,q,a0]=cvqm(fix(m),q,a0);
elseif(q<=50.0 & kd==2);
a0=(((-1.411114d-5.*q+9.730514d-4).*q-3.097887d-3).*q+3.533597d-2).*q+49.0547;
elseif(q<=40.0 & kd==3);
a0=((-3.043872d-4.*q+2.05511d-2).*q-9.16292d-2).*q+49.19035;
else;
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
end;
elseif(m>=8);
if(q<=3..*m);
[m,q,a0]=cvqm(fix(m),q,a0);
elseif(q>m.*m);
[kd,m,q,a0]=cvql(fix(kd),fix(m),q,a0);
elseif(m==8 & kd==1);
a0=(((8.634308d-6.*q-2.100289d-3).*q+.169072).*q-4.64336).*q+109.4211;
elseif(m==8 & kd==4);
a0=((-6.7842d-5.*q+2.2057d-3).*q+.48296).*q+56.59;
elseif(m==9 & kd==2);
a0=(((2.906435d-6.*q-1.019893d-3).*q+.1101965).*q-3.821851).*q+127.6098;
elseif(m==9 & kd==3);
a0=((-9.577289d-5.*q+.01043839).*q+.06588934).*q+78.0198;
elseif(m==10 & kd==1);
a0=(((5.44927d-7.*q-3.926119d-4).*q+.0612099).*q-2.600805).*q+138.1923;
elseif(m==10 & kd==4);
a0=((-7.660143d-5.*q+.01132506).*q-.09746023).*q+99.29494;
elseif(m==11 & kd==2);
a0=(((-5.67615d-7.*q+7.152722d-6).*q+.01920291).*q-1.081583).*q+140.88;
elseif(m==11 & kd==3);
a0=((-6.310551d-5.*q+.0119247).*q-.2681195).*q+123.667;
elseif(m==12 & kd==1);
a0=(((-2.38351d-7.*q-2.90139d-5).*q+.02023088).*q-1.289).*q+171.2723;
elseif(m==12 & kd==4);
a0=(((3.08902d-7.*q-1.577869d-4).*q+.0247911).*q-1.05454).*q+161.471;
end;
end;
end
%--/cvql.f90  processed by SPAG 6.53Rc at 13:11 on 28 Jun 2004
%
%
function [kd,m,q,a0]=cvql(kd,m,q,a0,varargin);
%
%     ========================================================
%     Purpose: Compute the characteristic value of Mathieu
%     functions  for q ò 3m
%     Input :  m  --- Order of Mathieu functions
%     q  --- Parameter of Mathieu functions
%     Output:  A0 --- Initial characteristic value
%     ========================================================
%
%
 c1=0;
 cv1=0;
 cv2=0;
 d1=0;
 d2=0;
 d3=0;
 d4=0;
 p1=0;
 p2=0;
 w=0;
 w2=0;
w3=0;
 w4=0;
 w6=0;
if(kd==1 | kd==2);
w=2.0d0.*fix(m)+1.0d0;
end;
if(kd==3 | kd==4);
w=2.0d0.*fix(m)-1.0d0;
end;
w2=w.*w;
w3=w.*w2;
w4=w2.*w2;
w6=w2.*w4;
d1=5.0+34.0./w2+9.0./w4;
d2=(33.0+410.0./w2+405.0./w4)./w;
d3=(63.0+1260.0./w2+2943.0./w4+486.0./w6)./w2;
d4=(527.0+15617.0./w2+69001.0./w4+41607.0./w6)./w3;
c1=128.0;
p2=q./w4;
p1=sqrt(p2);
cv1=-2.0.*q+2.0.*w.*sqrt(q)-(w2+1.0)./8.0;
cv2=(w+3.0./w)+d1./(32.0.*p1)+d2./(8.0.*c1.*p2);
cv2=cv2+d3./(64.0.*c1.*p1.*p2)+d4./(16.0.*c1.*c1.*p2.*p2);
a0=cv1-cv2./(c1.*p1);
end
%--/cvqm.f90  processed by SPAG 6.53Rc at 13:11 on 28 Jun 2004
%
%
function [m,q,a0]=cvqm(m,q,a0,varargin);
%
%     =====================================================
%     Purpose: Compute the characteristic value of Mathieu
%     functions for q ó m*m
%     Input :  m  --- Order of Mathieu functions
%     q  --- Parameter of Mathieu functions
%     Output:  A0 --- Initial characteristic value
%     =====================================================
%
%
 hm1=0;
 hm3=0;
 hm5=0;
hm1=.5.*q./(fix(m).*fix(m)-1.0);
hm3=.25.*hm1.^3./(fix(m).*fix(m)-4.0);
hm5=hm1.*hm3.*q./((fix(m).*fix(m)-1.0).*(fix(m).*fix(m)-9.0));
a0=fix(m).*fix(m)+q.*(hm1+(5.0.*fix(m).*fix(m)+7.0).*hm3+(9.0.*fix(m).^4+58.0.*fix(m).*fix(m)+29.0).*hm5);
end
%--/jynb.f90  processed by SPAG 6.53Rc at 13:11 on 28 Jun 2004
%
%
function [n,x,nm,bj,dj,by,dy]=jynb(n,x,nm,bj,dj,by,dy,varargin);
%
%     =====================================================
%     Purpose: Compute Bessel functions Jn(x), Yn(x)and
%     their derivatives
%     Input :  x --- Argument of Jn(x)and Yn(x,x ò 0)
%     n --- Order of Jn(x)and Yn(x)
%     Output:  BJ(n)--- Jn(x)
%     DJ(n)--- Jn'(x)
%     BY(n)--- Yn(x)
%     DY(n)--- Yn'(x)
%     NM --- Highest order computed
%     Routines called:
%     MSTA1 and MSTA2 to calculate the starting
%     point for backward recurrence
%     =====================================================
%
%
 a=zeros(1,4);
 a1=zeros(1,4);
 b=zeros(1,4);
 b1=zeros(1,4);
 bj0=0;
 bj1=0;
 bjk=0;
 bs=0;
by0=0;
 by1=0;
 byk=0;
 cu=0;
 ec=0;
 f=0;
 f1=0;
 f2=0;
 p0=0;
 p1=0;
pi=0;
 q0=0;
 q1=0;
 r2p=0;
 s0=0;
 su=0;
 sv=0;
 t1=0;
 t2=0;
 k=0;
 m=0;
pi=3.141592653589793d0;
r2p=.63661977236758d0;
nm=fix(fix(n));
if(x<1.0d-100);
for k=0 : n;
bj(k+1)=0.0d0;
dj(k+1)=0.0d0;
by(k+1)=-1.0d+300;
dy(k+1)=1.0d+300;
end; k=fix( fix(n)+1);
bj(0+1)=1.0d0;
dj(1+1)=0.5d0;
return;
end;
if(x<=300.0 | n>fix(0.9.*x));
if(n==0);
nm=1;
end;
m=fix(msta1(x,200));
if(m<nm);
nm=fix(fix(m));
else;
m=fix(msta2(x,fix(nm),15));
end;
bs=0.0d0;
su=0.0d0;
sv=0.0d0;
f2=0.0d0;
f1=1.0d-100;
for k=m : -1: 0 ;
f=2.0d0.*(fix(k)+1.0d0)./x.*f1-f2;
if(k<=nm);
bj(k+1)=f;
end;
if(k==2.*fix(k./2)& k~=0);
bs=bs+2.0d0.*f;
su=su+(-1).^(fix(k)./2).*f./fix(k);
elseif(k>1);
sv=sv+(-1).^(fix(fix(k)./2)).*fix(k)./(fix(k).*fix(k)-1.0).*f;
end;
f2=f1;
f1=f;
end; k=fix( 0 -1);
s0=bs+f;
for k=0 : nm;
bj(k+1)=bj(fix(k)+1)./s0;
end; k=fix( fix(nm)+1);
ec=log(x./2.0d0)+0.5772156649015329d0;
by0=r2p.*(ec.*bj(0+1)-4.0d0.*su./s0);
by(0+1)=by0;
by1=r2p.*((ec-1.0d0).*bj(1+1)-bj(0+1)./x-4.0d0.*sv./s0);
by(1+1)=by1;
else;
a(:)=[-.7031250000000000d-01,.1121520996093750d+00,-.5725014209747314d+00,.6074042001273483d+01];
b(:)=[.7324218750000000d-01,-.2271080017089844d+00,.1727727502584457d+01,-.2438052969955606d+02];
a1(:)=[.1171875000000000d+00,-.1441955566406250d+00,.6765925884246826d+00,-.6883914268109947d+01];
b1(:)=[-.1025390625000000d+00,.2775764465332031d+00,-.1993531733751297d+01,.2724882731126854d+02];
t1=x-0.25d0.*pi;
p0=1.0d0;
q0=-0.125d0./x;
for k=1 : 4;
p0=p0+a(fix(k)).*x.^(-2.*fix(k));
q0=q0+b(fix(k)).*x.^(-2.*fix(k)-1);
end; k=fix( 4+1);
cu=sqrt(r2p./x);
bj0=cu.*(p0.*cos(t1)-q0.*sin(t1));
by0=cu.*(p0.*sin(t1)+q0.*cos(t1));
bj(0+1)=bj0;
by(0+1)=by0;
t2=x-0.75d0.*pi;
p1=1.0d0;
q1=0.375d0./x;
for k=1 : 4;
p1=p1+a1(fix(k)).*x.^(-2.*fix(k));
q1=q1+b1(fix(k)).*x.^(-2.*fix(k)-1);
end; k=fix( 4+1);
bj1=cu.*(p1.*cos(t2)-q1.*sin(t2));
by1=cu.*(p1.*sin(t2)+q1.*cos(t2));
bj(1+1)=bj1;
by(1+1)=by1;
for k=2 : nm;
bjk=2.0d0.*(fix(k)-1.0d0)./x.*bj1-bj0;
bj(k+1)=bjk;
bj0=bj1;
bj1=bjk;
end; k=fix( fix(nm)+1);
end;
dj(0+1)=-bj(1+1);
for k=1 : nm;
dj(k+1)=bj(fix(k)-1+1)-fix(k)./x.*bj(fix(k)+1);
end; k=fix( fix(nm)+1);
for k=2 : nm;
byk=2.0d0.*(fix(k)-1.0d0).*by1./x-by0;
by(k+1)=byk;
by0=by1;
by1=byk;
end; k=fix( fix(nm)+1);
dy(0+1)=-by(1+1);
for k=1 : nm;
dy(k+1)=by(fix(k)-1+1)-fix(k).*by(fix(k)+1)./x;
end; k=fix( fix(nm)+1);
end
%--/msta1.f90  processed by SPAG 6.53Rc at 13:11 on 28 Jun 2004
%
%
function [msta1Result]=msta1(x,mp,varargin);
%
%     ===================================================
%     Purpose: Determine the starting point for backward
%     recurrence such that the magnitude of
%     Jn(x)at that point is about 10^(-MP)
%     Input :  x     --- Argument of Jn(x)
%     MP    --- Value of magnitude
%     Output:  MSTA1 --- Starting point
%     ===================================================
%
%
 a0=0;
 f=0;
 f0=0;
 f1=0;
 it=0;
 n0=0;
 n1=0;
 nn=0;
a0=abs(x);
n0=fix(fix(1.1.*a0)+1);
f0=envj(fix(n0),a0)-fix(mp);
n1=fix(fix(n0)+5);
f1=envj(fix(n1),a0)-fix(mp);
for it=1 : 20;
nn=fix(fix(n1)-(fix(n1)-fix(n0))./(1.0d0-f0./f1));
f=envj(fix(nn),a0)-fix(mp);
if(abs(nn-n1)<1);
break;
end;
n0=fix(fix(n1));
f0=f1;
n1=fix(fix(nn));
f1=f;
end;
msta1Result=fix(fix(nn));
end
%--/msta2.f90  processed by SPAG 6.53Rc at 13:11 on 28 Jun 2004
%
%
function [msta2Result]=msta2(x,n,mp,varargin);
%
%     ===================================================
%     Purpose: Determine the starting point for backward
%     recurrence such that all Jn(x)has MP
%     significant digits
%     Input :  x  --- Argument of Jn(x)
%     n  --- Order of Jn(x)
%     MP --- Significant digit
%     Output:  MSTA2 --- Starting point
%     ===================================================
%
%
 a0=0;
 ejn=0;
 f=0;
 f0=0;
 f1=0;
 hmp=0;
 obj=0;
 it=0;
 n0=0;
 n1=0;
 nn=0;
a0=abs(x);
hmp=0.5d0.*fix(mp);
ejn=envj(fix(n),a0);
if(ejn<=hmp);
obj=fix(mp);
n0=fix(fix(1.1.*a0));
else;
obj=hmp+ejn;
n0=fix(fix(n));
end;
f0=envj(fix(n0),a0)-obj;
n1=fix(fix(n0)+5);
f1=envj(fix(n1),a0)-obj;
for it=1 : 20;
nn=fix(fix(n1)-(fix(n1)-fix(n0))./(1.0d0-f0./f1));
f=envj(fix(nn),a0)-obj;
if(abs(nn-n1)<1);
break;
end;
n0=fix(fix(n1));
f0=f1;
n1=fix(fix(nn));
f1=f;
end;
msta2Result=fix(fix(nn)+10);
end
%--/envj.f90  processed by SPAG 6.53Rc at 13:11 on 28 Jun 2004
%
function [envjResult]=envj(n,x,varargin);
%
envjResult=0.5d0.*log10(6.28d0.*fix(n))-fix(n).*log10(1.36d0.*x./fix(n));
end

