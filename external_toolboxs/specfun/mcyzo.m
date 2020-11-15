function mcyzo
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     ===========================================================
%     Purpose : This program evaluates the complex zeros of
%     Y0(z), Y0'(z), Y1(z)and Y1'(z), and their
%     associated values at the zeros using the
%     modified Newton's iteration method
%     Input:    NT --- Total number of roots/zeros
%     KF --- Function choice code
%     KF=0 for  Y0(z)& Y1(z0)
%     KF=1 for  Y1(z)& Y0(z1)
%     KF=2 for  Y1'(z)& Y1(z1')
%     KC --- Choice code
%     KC=0 for complex roots
%     KC=1 for real roots
%     Output:   ZO(L)--- L-th zero of Y0(z)or Y1(z)or Y1'(z)
%     ZV(L)--- Value of Y0'(z)or Y1'(z)or Y1(z)
%     at the L-th zero
%     Examples: NT = 5
%     No.      z0, Zeros of Y0(z)Y1(z0)
%     -----------------------------------------------------------------
%     1   -2.403016632 + i .5398823130   .1007476893 - i .8819677101
%     2   -5.519876702 + i .5471800106  -.0292464182 + i .5871695027
%     3   -8.653672403 + i .5484120673   .0149080637 - i .4694587524
%     4  -11.791512030 + i .5488191184  -.0093736817 + i .4023045429
%     5  -14.930906564 + i .5490008289   .0065788031 - i .3575673214
%     No.      z1, Zeros of Y1(z)Y0(z1)
%     -----------------------------------------------------------------
%     1    -.502743273 + i .7862437145  -.4595276847 + i1.3171019361
%     2   -3.833535193 + i .5623565382   .0483019087 - i .6925128842
%     3   -7.015903683 + i .5533930459  -.0201269494 + i .5186425332
%     4  -10.173573834 + i .5512733877   .0116140017 - i .4320329636
%     5  -13.323739307 + i .5504585830  -.0077719300 + i .3779698048
%     No.      z1', Zeros of Y1'(z)Y1(z1')
%     ----------------------------------------------------------------
%     1     .576785129 + i .9039847922  -.7634970879 + i .5892448647
%     2   -1.940477342 - i .7211859189   .1620640057 + i .9520278864
%     3   -5.333478617 - i .5672196368  -.0317940081 - i .5968536736
%     4   -8.536768577 - i .5560607040   .0154177166 + i .4726011652
%     5  -11.706175219 - i .5528590607  -.0095443768 - i .4037533396
%     ============================================================
nt=[];kf=[];kc=[];zo=[];zv=[];
 zo=zeros(1,50);
zv=zeros(1,50);
fprintf(1,'%s \n','please enter nt, kf and kc');
fprintf(1,'%s \n','  nt --- total number of the roots');
fprintf(1,'%s \n','  kf  --- function choice code');
fprintf(1,'%s \n','          kf=0 for y0(z)& y1(z0)');
fprintf(1,'%s \n','          kf=1 for y1(z)& y0(z1)');
fprintf(1,'%s \n','          kf=2 for y1''(z)& y1(z1'')');
fprintf(1,'%s \n','  kc  --- choice code');
fprintf(1,'%s \n','          kc=0 for complex roots');
fprintf(1,'%s \n','          kc=1 for real roots');
%     READ(*,*)NT,KF,KC
nt=5;
kf=1;
kc=0;
fprintf(1,[repmat(' ',1,2),'nt=','%3g',',  ','kf=','%3g',',  ','kc=','%3g' ' \n'],nt,kf,kc);
fprintf(1,'%0.15g \n');
fprintf(1,[repmat(' ',1,20),'*****    please wait     *****' ' \n']);
[nt,kf,kc,zo,zv]=cyzo(nt,kf,kc,zo,zv);
fprintf(1,'%0.15g \n');
if(kf == 0);
fprintf(1,'%s ',' no.          z0, zeros of y0(z)');fprintf(1,'%s \n', '                 y1(z0)');
elseif(kf == 1);
fprintf(1,'%s ',' no.          z1, zeros of y1(z)');fprintf(1,'%s \n', '                 y0(z1)');
elseif(kf == 2);
fprintf(1,'%s ',' no.        z1'', zeros of y1''(z)');fprintf(1,'%s \n', '                y1(z1'')');
end;
fprintf(1,'%s ','--------------------------------------');fprintf(1,'%s \n','----------------------------');
for  i=1:nt;
fprintf(1,[repmat(' ',1,1),'%3g',repmat(' ',1,2),'%15.9g',repmat('%15.10g',1,3) ' \n'],i,zo(i),zv(i));
end;  i=nt+1;
%format(20x,'*****    please wait     *****');
%format(2x,',i3,',  ',',i3,',  ',',i3);
%format(1x,i3,2x,f15.9,3f15.10);
end
function [nt,kf,kc,zo,zv]=cyzo(nt,kf,kc,zo,zv,varargin);
%     ===========================================================
%     Purpose : Compute the complex zeros of Y0(z), Y1(z)and
%     Y1'(z), and their associated values at the zeros
%     using the modified Newton's iteration method
%     Input:    NT --- Total number of zeros/roots
%     KF --- Function choice code
%     KF=0 for  Y0(z)& Y1(z0)
%     KF=1 for  Y1(z)& Y0(z1)
%     KF=2 for  Y1'(z)& Y1(z1')
%     KC --- Choice code
%     KC=0 for complex roots
%     KC=1 for real roots
%     Output:   ZO(L)--- L-th zero of Y0(z)or Y1(z)or Y1'(z)
%     ZV(L)--- Value of Y0'(z)or Y1'(z)or Y1(z)
%     at the L-th zero
%     Routine called: CY01 for computing Y0(z)and Y1(z), and
%     their derivatives
%     ===========================================================
z=[];zf=[];zd=[];
w=0.0;
if(kc == 0);
x=-2.4d0;
y=0.54d0;
h=3.14;
elseif(kc == 1);
x=0.89;
y=0.0;
h=-3.14;
end;
if(kf == 1)x=-0.503; end;
if(kf == 2)x=0.577; end;
zero=complex(x,y);
z=zero;
for  nr=1:nt;
if(nr ~= 1)z=zo(nr-1)-h; end;
it=0;
while (1);
it=it+1;
[kf,z,zf,zd]=cy01(fix(kf),z,zf,zd);
zp=complex(1.0d0,0.0d0);
for  i=1:nr-1;
zp=zp.*(z-zo(i));
end;  i=nr-1+1;
zfd=zf./zp;
zq=complex(0.0d0,0.0d0);
for  i=1:nr-1;
zw=complex(1.0d0,0.0d0);
for  j=1:nr-1;
if(~(j == i));
zw=zw.*(z-zo(j));
end;
end;  j=nr-1+1;
zq=zq+zw;
end;  i=nr-1+1;
zgd=(zd-zq.*zfd)./zp;
z=z-zfd./zgd;
w0=w;
w=abs(z);
if(~(it <= 50&abs((w-w0)./w)> 1.0d-12))break; end;
end;
zo(nr)=z;
end;
for  i=1:nt;
z=zo(i);
if(kf == 0|kf == 2);
[dumvar1,z,zf,zd]=cy01(1,z,zf,zd);
zv(i)=zf;
elseif(kf == 1);
[dumvar1,z,zf,zd]=cy01(0,z,zf,zd);
zv(i)=zf;
end;
end;  i=fix(nt)+1;
return;
end
function [kf,z,zf,zd]=cy01(kf,z,zf,zd,varargin);
%     ===========================================================
%     Purpose: Compute complex Bessel functions Y0(z), Y1(z)
%     and their derivatives
%     Input :  z  --- Complex argument of Yn(z,n=0,1)
%     KF --- Function choice code
%     KF=0 for ZF=Y0(z)and ZD=Y0'(z)
%     KF=1 for ZF=Y1(z)and ZD=Y1'(z)
%     KF=2 for ZF=Y1'(z)and ZD=Y1''(z)
%     Output:  ZF --- Y0(z)or Y1(z)or Y1'(z)
%     ZD --- Y0'(z)or Y1'(z)or Y1''(z)
%     ===========================================================
 a=zeros(1,12);
b=zeros(1,12);
a1=zeros(1,12);
b1=zeros(1,12);
pi=3.141592653589793d0;
el=0.5772156649015329d0;
rp2=2.0d0./pi;
ci=complex(0.0d0,1.0d0);
a0=abs(z);
z2=z.*z;
z1=z;
if(a0 == 0.0d0);
cbj0=complex(1.0d0,0.0d0);
cbj1=complex(0.0d0,0.0d0);
cby0=-complex(1.0d300,0.0d0);
cby1=-complex(1.0d300,0.0d0);
cdy0=complex(1.0d300,0.0d0);
cdy1=complex(1.0d300,0.0d0);
else;
if(real(z)< 0.0)z1=-z; end;
if(a0 <= 12.0);
cbj0=complex(1.0d0,0.0d0);
cr=complex(1.0d0,0.0d0);
for  k=1:40;
cr=-0.25d0.*cr.*z2./(k.*k);
cbj0=cbj0+cr;
if(abs(cr)< abs(cbj0).*1.0d-15)break; end;
end;
cbj1=complex(1.0d0,0.0d0);
cr=complex(1.0d0,0.0d0);
for  k=1:40;
cr=-0.25d0.*cr.*z2./(k.*(k+1.0d0));
cbj1=cbj1+cr;
if(abs(cr)< abs(cbj1).*1.0d-15)break; end;
end;
cbj1=0.5d0.*z1.*cbj1;
w0=0.0d0;
cr=complex(1.0d0,0.0d0);
cs=complex(0.0d0,0.0d0);
for  k=1:40;
w0=w0+1.0d0./k;
cr=-0.25d0.*cr./(k.*k).*z2;
cp=cr.*w0;
cs=cs+cp;
if(abs(cp)< abs(cs).*1.0d-15)break; end;
end;
cby0=rp2.*(log(z1./2.0d0)+el).*cbj0-rp2.*cs;
w1=0.0d0;
cr=complex(1.0d0,0.0d0);
cs=complex(1.0d0,0.0d0);
for  k=1:40;
w1=w1+1.0d0./k;
cr=-0.25d0.*cr./(k.*(k+1)).*z2;
cp=cr.*(2.0d0.*w1+1.0d0./(k+1.0d0));
cs=cs+cp;
if(abs(cp)< abs(cs).*1.0d-15)break; end;
end;
cby1=rp2.*((log(z1./2.0d0)+el).*cbj1-1.0d0./z1-.25d0.*z1.*cs);
else;
a(:)=[-.703125d-01,.112152099609375d+00,-.5725014209747314d+00,.6074042001273483d+01,-.1100171402692467d+03,.3038090510922384d+04,-.1188384262567832d+06,.6252951493434797d+07,-.4259392165047669d+09,.3646840080706556d+11,-.3833534661393944d+13,.4854014686852901d+15];
b(:)=[.732421875d-01,-.2271080017089844d+00,.1727727502584457d+01,-.2438052969955606d+02,.5513358961220206d+03,-.1825775547429318d+05,.8328593040162893d+06,-.5006958953198893d+08,.3836255180230433d+10,-.3649010818849833d+12,.4218971570284096d+14,-.5827244631566907d+16];
a1(:)=[.1171875d+00,-.144195556640625d+00,.6765925884246826d+00,-.6883914268109947d+01,.1215978918765359d+03,-.3302272294480852d+04,.1276412726461746d+06,-.6656367718817688d+07,.4502786003050393d+09,-.3833857520742790d+11,.4011838599133198d+13,-.5060568503314727d+15];
b1(:)=[-.1025390625d+00,.2775764465332031d+00,-.1993531733751297d+01,.2724882731126854d+02,-.6038440767050702d+03,.1971837591223663d+05,-.8902978767070678d+06,.5310411010968522d+08,-.4043620325107754d+10,.3827011346598605d+12,-.4406481417852278d+14,.6065091351222699d+16];
k0=12;
if(a0 >= 35.0)k0=10; end;
if(a0 >= 50.0)k0=8; end;
ct1=z1-.25d0.*pi;
cp0=complex(1.0d0,0.0d0);
for  k=1:k0;
cp0=cp0+a(k).*z1.^(-2.*k);
end;  k=k0+1;
cq0=-0.125d0./z1;
for  k=1:k0;
cq0=cq0+b(k).*z1.^(-2.*k-1);
end;  k=k0+1;
cu=sqrt(rp2./z1);
cbj0=cu.*(cp0.*cos(ct1)-cq0.*sin(ct1));
cby0=cu.*(cp0.*sin(ct1)+cq0.*cos(ct1));
ct2=z1-.75d0.*pi;
cp1=complex(1.0d0,0.0d0);
for  k=1:k0;
cp1=cp1+a1(k).*z1.^(-2.*k);
end;  k=k0+1;
cq1=0.375d0./z1;
for  k=1:k0;
cq1=cq1+b1(k).*z1.^(-2.*k-1);
end;  k=k0+1;
cbj1=cu.*(cp1.*cos(ct2)-cq1.*sin(ct2));
cby1=cu.*(cp1.*sin(ct2)+cq1.*cos(ct2));
end;
if(real(z)< 0.0);
if(imag(z)< 0.0)cby0=cby0-2.0d0.*ci.*cbj0; end;
if(imag(z)> 0.0)cby0=cby0+2.0d0.*ci.*cbj0; end;
if(imag(z)< 0.0)cby1=-(cby1-2.0d0.*ci.*cbj1); end;
if(imag(z)> 0.0)cby1=-(cby1+2.0d0.*ci.*cbj1); end;
cbj1=-cbj1;
end;
cdy0=-cby1;
cdy1=cby0-1.0d0./z.*cby1;
end;
if(kf == 0);
zf=cby0;
zd=cdy0;
elseif(kf == 1);
zf=cby1;
zd=cdy1;
elseif(kf == 2);
zf=cdy1;
zd=-cdy1./z-(1.0d0-1.0d0./(z.*z)).*cby1;
end;
return;
end

