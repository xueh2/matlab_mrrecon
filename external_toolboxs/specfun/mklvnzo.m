function mklvnzo
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     ==============================================================
%     Purpose: This program computes the first NT zeros of Kelvin
%     functions and their derivatives using subroutine
%     KLVNZO
%     Input :  NT --- Total number of zeros
%     Example: NT = 5
%     Zeros of Kelvin functions ber x, bei x, ker x and kei x
%     m       ber x          bei x          ker x          kei x
%     ---------------------------------------------------------------
%     1     2.84891782     5.02622395     1.71854296     3.91466761
%     2     7.23882945     9.45540630     6.12727913     8.34422506
%     3    11.67396355    13.89348785    10.56294271    12.78255715
%     4    16.11356383    18.33398346    15.00268812    17.22314372
%     5    20.55463158    22.77543929    19.44381663    21.66464214
%     Zeros of Kelvin Functions ber'x, bei'x, ker'x and kei'x
%     m       ber'x          bei'x          ker'x          kei'x
%     ---------------------------------------------------------------
%     1     6.03871081     3.77267330     2.66583979     4.93181194
%     2    10.51364251     8.28098785     7.17212212     9.40405458
%     3    14.96844542    12.74214752    11.63218639    13.85826916
%     4    19.41757493    17.19343175    16.08312025    18.30717294
%     5    23.86430432    21.64114394    20.53067845    22.75379258
%     ==============================================================
nt=[];r1=[];r2=[];r3=[];r4=[];r5=[];r6=[];r7=[];r8=[];
 r1=zeros(1,50);
r2=zeros(1,50);
r3=zeros(1,50);
r4=zeros(1,50);
r5=zeros(1,50);
r6=zeros(1,50);
 r7=zeros(1,50);
r8=zeros(1,50);
fprintf(1,[repmat(' ',1,1),'nt is the number of the zeros' ' \n']);
fprintf(1,'%s \n','please enter nt ');
%     READ(*,*)NT
nt=5;
fprintf(1,[repmat(' ',1,4),'zeros of kelvin functions ber x, bei x,',' ker x and kei x' ' \n']);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   m       ber x          bei x          ker x');fprintf(1,'%s \n', '          kei x');
fprintf(1,'%s ',' ------------------------------------------------');fprintf(1,'%s \n','---------------');
[nt,dumvar2,r1]=klvnzo(nt,1,r1);
[nt,dumvar2,r2]=klvnzo(nt,2,r2);
[nt,dumvar2,r3]=klvnzo(nt,3,r3);
[nt,dumvar2,r4]=klvnzo(nt,4,r4);
for  l=1:nt;
fprintf(1,[repmat(' ',1,1),'%3g',repmat(' ',1,1),'%14.8g',repmat(' ',1,1),'%14.8g',repmat(' ',1,1),'%14.8g',repmat(' ',1,1),'%14.8g' ' \n'],l,r1(l),r2(l),r3(l),r4(l));
end;  l=nt+1;
[nt,dumvar2,r5]=klvnzo(nt,5,r5);
[nt,dumvar2,r6]=klvnzo(nt,6,r6);
[nt,dumvar2,r7]=klvnzo(nt,7,r7);
[nt,dumvar2,r8]=klvnzo(nt,8,r8);
fprintf(1,'%0.15g \n');
fprintf(1,[repmat(' ',1,4),'zeros of kelvin functions ber''x, bei''x,',' ker''x and kei''x' ' \n']);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','   m       ber''x          bei''x          ker''x');fprintf(1,'%s \n','          kei''x');
fprintf(1,'%s ',' ------------------------------------------------');fprintf(1,'%s \n','---------------');
for  l=1:nt;
fprintf(1,[repmat(' ',1,1),'%3g',repmat(' ',1,1),'%14.8g',repmat(' ',1,1),'%14.8g',repmat(' ',1,1),'%14.8g',repmat(' ',1,1),'%14.8g' ' \n'],l,r5(l),r6(l),r7(l),r8(l));
end;  l=nt+1;
%format(1x,i3,1x,f14.8,1x,f14.8,1x,f14.8,1x,f14.8);
%format(4x,'zeros of kelvin functions ber x, bei x,' ,' ker x and kei x');
%format(4x,'zeros of kelvin functions ber''x, bei''x,',' ker''x and kei''x');
%format(1x,'nt is the number of the zeros');
end
function [nt,kd,zo]=klvnzo(nt,kd,zo,varargin);
%     ====================================================
%     Purpose: Compute the zeros of Kelvin functions
%     Input :  NT  --- Total number of zeros
%     KD  --- Function code
%     KD=1 to 8 for ber x, bei x, ker x, kei x,
%     ber'x, bei'x, ker'x and kei'x,
%     respectively.
%     Output:  ZO(M)--- the M-th zero of Kelvin function
%     for code KD
%     Routine called:
%     KLVNA for computing Kelvin functions and
%     their derivatives
%     ====================================================
rt=[];ber=[];bei=[];ger=[];gei=[];der=[];dei=[];her=[];hei=[];
  rt0=zeros(1,8);
rt0(1)=2.84891;
rt0(2)=5.02622;
rt0(3)=1.71854;
rt0(4)=3.91467;
rt0(5)=6.03871;
rt0(6)=3.77268;
rt0(7)=2.66584;
rt0(8)=4.93181;
rt=rt0(fix(kd));
for  m=1:nt;
while (1);
[rt,ber,bei,ger,gei,der,dei,her,hei]=klvna(rt,ber,bei,ger,gei,der,dei,her,hei);
if(kd == 1);
rt=rt-ber./der;
elseif(kd == 2);
rt=rt-bei./dei;
elseif(kd == 3);
rt=rt-ger./her;
elseif(kd == 4);
rt=rt-gei./hei;
elseif(kd == 5);
ddr=-bei-der./rt;
rt=rt-der./ddr;
elseif(kd == 6);
ddi=ber-dei./rt;
rt=rt-dei./ddi;
elseif(kd == 7);
gdr=-gei-her./rt;
rt=rt-her./gdr;
else;
gdi=ger-hei./rt;
rt=rt-hei./gdi;
end;
if(abs(rt-rt0(kd))> 5.0d-10);
rt0(kd)=rt;
else;
break;
end;
end;
zo(m)=rt;
rt=rt+4.44d0;
end;
return;
end
function [x,ber,bei,ger,gei,der,dei,her,hei]=klvna(x,ber,bei,ger,gei,der,dei,her,hei,varargin);
%     ======================================================
%     Purpose: Compute Kelvin functions ber x, bei x, ker x
%     and kei x, and their derivatives(x > 0)
%     Input :  x   --- Argument of Kelvin functions
%     Output:  BER --- ber x
%     BEI --- bei x
%     GER --- ker x
%     GEI --- kei x
%     DER --- ber'x
%     DEI --- bei'x
%     HER --- ker'x
%     HEI --- kei'x
%     ================================================
pi=3.141592653589793d0;
el=.5772156649015329d0;
eps=1.0d-15;
if(x == 0.0d0);
ber=1.0d0;
bei=0.0d0;
ger=1.0d+300;
gei=-.25d0.*pi;
der=0.0d0;
dei=0.0d0;
her=-1.0d+300;
hei=0.0d0;
return;
end;
x2=.25d0.*x.*x;
x4=x2.*x2;
if(abs(x)< 10.0d0);
ber=1.0d0;
r=1.0d0;
for  m=1:60;
r=-.25d0.*r./(m.*m)./(2.0d0.*m-1.0d0).^2.*x4;
ber=ber+r;
if(abs(r./ber)< eps)break; end;
end;
bei=x2;
r=x2;
for  m=1:60;
r=-.25d0.*r./(m.*m)./(2.0d0.*m+1.0d0).^2.*x4;
bei=bei+r;
if(abs(r./bei)< eps)break; end;
end;
ger=-(log(x./2.0d0)+el).*ber+.25d0.*pi.*bei;
r=1.0d0;
gs=0.0d0;
for  m=1:60;
r=-.25d0.*r./(m.*m)./(2.0d0.*m-1.0d0).^2.*x4;
gs=gs+1.0d0./(2.0d0.*m-1.0d0)+1.0d0./(2.0d0.*m);
ger=ger+r.*gs;
if(abs(r.*gs./ger)< eps)break; end;
end;
gei=x2-(log(x./2.0d0)+el).*bei-.25d0.*pi.*ber;
r=x2;
gs=1.0d0;
for  m=1:60;
r=-.25d0.*r./(m.*m)./(2.0d0.*m+1.0d0).^2.*x4;
gs=gs+1.0d0./(2.0d0.*m)+1.0d0./(2.0d0.*m+1.0d0);
gei=gei+r.*gs;
if(abs(r.*gs./gei)< eps)break; end;
end;
der=-.25d0.*x.*x2;
r=der;
for  m=1:60;
r=-.25d0.*r./m./(m+1.0d0)./(2.0d0.*m+1.0d0).^2.*x4;
der=der+r;
if(abs(r./der)< eps)break; end;
end;
dei=.5d0.*x;
r=dei;
for  m=1:60;
r=-.25d0.*r./(m.*m)./(2.d0.*m-1.d0)./(2.d0.*m+1.d0).*x4;
dei=dei+r;
if(abs(r./dei)< eps)break; end;
end;
r=-.25d0.*x.*x2;
gs=1.5d0;
her=1.5d0.*r-ber./x-(log(x./2.d0)+el).*der+.25.*pi.*dei;
for  m=1:60;
r=-.25d0.*r./m./(m+1.0d0)./(2.0d0.*m+1.0d0).^2.*x4;
gs=gs+1.0d0./(2.*m+1.0d0)+1.0d0./(2.*m+2.0d0);
her=her+r.*gs;
if(abs(r.*gs./her)< eps)break; end;
end;
r=.5d0.*x;
gs=1.0d0;
hei=.5d0.*x-bei./x-(log(x./2.d0)+el).*dei-.25.*pi.*der;
for  m=1:60;
r=-.25d0.*r./(m.*m)./(2.*m-1.0d0)./(2.*m+1.0d0).*x4;
gs=gs+1.0d0./(2.0d0.*m)+1.0d0./(2.*m+1.0d0);
hei=hei+r.*gs;
if(abs(r.*gs./hei)< eps)return; end;
end;  m=60+1;
else;
pp0=1.0d0;
pn0=1.0d0;
qp0=0.0d0;
qn0=0.0d0;
r0=1.0d0;
km=18;
if(abs(x)>= 40.0)km=10; end;
fac=1.0d0;
for  k=1:km;
fac=-fac;
xt=.25d0.*k.*pi-fix(.125d0.*k).*2.0d0.*pi;
cs=cos(xt);
ss=sin(xt);
r0=.125d0.*r0.*(2.0d0.*k-1.0d0).^2./k./x;
rc=r0.*cs;
rs=r0.*ss;
pp0=pp0+rc;
pn0=pn0+fac.*rc;
qp0=qp0+rs;
qn0=qn0+fac.*rs;
end;  k=km+1;
xd=x./sqrt(2.0d0);
xe1=exp(xd);
xe2=exp(-xd);
xc1=1.d0./sqrt(2.0d0.*pi.*x);
xc2=sqrt(.5d0.*pi./x);
cp0=cos(xd+.125d0.*pi);
cn0=cos(xd-.125d0.*pi);
sp0=sin(xd+.125d0.*pi);
sn0=sin(xd-.125d0.*pi);
ger=xc2.*xe2.*(pn0.*cp0-qn0.*sp0);
gei=xc2.*xe2.*(-pn0.*sp0-qn0.*cp0);
ber=xc1.*xe1.*(pp0.*cn0+qp0.*sn0)-gei./pi;
bei=xc1.*xe1.*(pp0.*sn0-qp0.*cn0)+ger./pi;
pp1=1.0d0;
pn1=1.0d0;
qp1=0.0d0;
qn1=0.0d0;
r1=1.0d0;
fac=1.0d0;
for  k=1:km;
fac=-fac;
xt=.25d0.*k.*pi-fix(.125d0.*k).*2.0d0.*pi;
cs=cos(xt);
ss=sin(xt);
r1=.125d0.*r1.*(4.d0-(2.0d0.*k-1.0d0).^2)./k./x;
rc=r1.*cs;
rs=r1.*ss;
pp1=pp1+fac.*rc;
pn1=pn1+rc;
qp1=qp1+fac.*rs;
qn1=qn1+rs;
end;  k=km+1;
her=xc2.*xe2.*(-pn1.*cn0+qn1.*sn0);
hei=xc2.*xe2.*(pn1.*sn0+qn1.*cn0);
der=xc1.*xe1.*(pp1.*cp0+qp1.*sp0)-hei./pi;
dei=xc1.*xe1.*(pp1.*sp0-qp1.*cp0)+her./pi;
end;
return;
end

