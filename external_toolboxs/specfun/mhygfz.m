function mhygfz
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%     ============================================================
%     Purpose: This program computes hypergeometric function for
%     a complex argument, F(a,b,c,z), using subroutine
%     HYGFZ
%     Input :  a --- Parameter
%     b --- Parameter
%     c --- Parameter,  c <> 0,-1,-2,...
%     x --- Real part of complex argument z
%     y --- Imaginary part of complex argument z
%(z = x+iy)
%     Output:  ZHF --- F(a,b,c,z)
%     Examples:
%     a     b     c    z = x+ iy             F(a,b,c,z)
%     --------------------------------------------------------------
%     3.2   1.8   6.7   1.0+0.0 i    .54689992D+01+.00000000D+00 i
%     3.2  -1.8   6.7   1.0+0.0 i    .33750635D+00+.00000000D+00 i
%     -5.0   3.3   6.7   5.2+4.8 i    .11682745D+03+.60389104D+03 i
%     3.3  -6.0   3.7   5.2-4.8 i    .17620425D+05+.38293812D+05 i
%     -7.0   3.3  -3.7   5.2-4.8 i   -.11772779D+11-.14382286D+11 i
%     4.3  -8.0  -3.7   5.2+4.8 i    .13161188D+13-.10129870D+12 i
%     3.3   5.8   6.7   0.2+0.1 i    .17330557D+01+.63401030D+00 i
%     3.5  -2.4   6.7   0.2+0.5 i    .64762241D+00-.52110507D+00 i
%     3.3   4.3   6.7   0.8+0.3 i   -.14830086D+01+.83744258D+01 i
%     7.0   5.0   4.1   3.0-1.0 i   -.40376095D-02-.29566326D-02 i
%     5.0   7.0   4.1   3.0-1.0 i   -.40376095D-02-.29566326D-02 i
%     3.5   1.2   9.7   0.6+0.9 i    .10343044D+01+.54473814D+00 i
%     2.1   5.4   9.7   0.5+0.7 i    .68850442D+00+.12274187D+01 i
%     8.7   3.2   6.7   0.5+0.7 i   -.90046505D+00-.11198900D+01 i
%     8.7   2.7   6.7   0.6+0.9 i   -.46083890D+00-.54575701D+00 i
%     ============================================================
a=[];b=[];c=[];z=[];zhf=[];
fprintf(1,'%s \n','please enter a,b,c,x and y ');
%     READ(*,*)A,B,C,X,Y
a=8.7;
b=2.7;
c=6.7;
x=.6;
y=.9;
z=complex(x,y);
[a,b,c,z,zhf]=hygfz(a,b,c,z,zhf);
fprintf(1,'%0.15g \n');
fprintf(1,'%s ','     a      b      c      x      y');fprintf(1,'%s \n','          re[f]im[f]');
fprintf(1,'%s ','   --------------------------------');fprintf(1,'%s \n','-----------------------------------');
fprintf(1,[repmat(' ',1,1),repmat('%7.1g',1,5),repmat(' ',1,2),repmat('%16.8g',1,2) ' \n'],a,b,c,x,y,zhf);
%format(1x,5f7.1,2x,2d16.8);
end
function [a,b,c,z,zhf]=hygfz(a,b,c,z,zhf,varargin);
%     ======================================================
%     Purpose: Compute the hypergeometric function for a
%     complex argument, F(a,b,c,z)
%     Input :  a --- Parameter
%     b --- Parameter
%     c --- Parameter,  c <> 0,-1,-2,...
%     z --- Complex argument
%     Output:  ZHF --- F(a,b,c,z)
%     Routines called:
%(1)GAMMA for computing gamma function
%(2)PSI for computing psi function
%     ======================================================
gc=[];gcab=[];gca=[];gcb=[];g1=[];g2=[];g3=[];ga=[];gb=[];m=[];gam=[];gbm=[];pa=[];pb=[];gabc=[];gab=[];gba=[];pca=[];pac=[];gm=[];t0=[];g0=[];k=[];gcbk=[];
  l0=false;
 l1=false;
 l2=false;
 l3=false;
 l4=false;
 l5=false;
 l6=false;
zw=0.0;
x=real(z);
y=imag(z);
eps=1.0d-15;
l0=c == fix(c)&c < 0.0d0;
l1=abs(1.0d0-x)< eps&y == 0.0d0&c-a-b <= 0.0d0;
l2=abs(z+1.0d0)< eps&abs(c-a+b-1.0d0)< eps;
l3=a == fix(a)&a < 0.0d0;
l4=b == fix(b)&b < 0.0d0;
l5=c-a == fix(c-a)&c-a <= 0.0d0;
l6=c-b == fix(c-b)&c-b <= 0.0d0;
aa=a;
bb=b;
a0=abs(z);
if(a0 > 0.95d0)eps=1.0d-8; end;
pi=3.141592653589793d0;
el=.5772156649015329d0;
if(l0|l1);
fprintf(1,'%s \n','the hypergeometric series is divergent');
return;
end;
if(a0 == 0.0d0|a == 0.0d0|b == 0.0d0);
zhf=complex(1.0d0,0.0d0);
elseif(z == 1.0d0&c-a-b > 0.0d0);
[c,gc]=gamma(c,gc);
[dumvar1,gcab]=gamma(c-a-b,gcab);
[dumvar1,gca]=gamma(c-a,gca);
[dumvar1,gcb]=gamma(c-b,gcb);
zhf=gc.*gcab./(gca.*gcb);
elseif(l2);
g0=sqrt(pi).*2.0d0.^(-a);
[c,g1]=gamma(c,g1);
[dumvar1,g2]=gamma(1.0d0+a./2.0d0-b,g2);
[dumvar1,g3]=gamma(0.5d0+0.5d0.*a,g3);
zhf=g0.*g1./(g2.*g3);
elseif(l3|l4);
if(l3)nm=fix(abs(a)); end;
if(l4)nm=fix(abs(b)); end;
zhf=complex(1.0d0,0.0d0);
zr=complex(1.0d0,0.0d0);
for  k=1:nm;
zr=zr.*(a+k-1.0d0).*(b+k-1.0d0)./(k.*(c+k-1.0d0)).*z;
zhf=zhf+zr;
end;  k=nm+1;
elseif(l5|l6);
if(l5)nm=fix(abs(c-a)); end;
if(l6)nm=fix(abs(c-b)); end;
zhf=complex(1.0d0,0.0d0);
zr=complex(1.0d0,0.0d0);
for  k=1:nm;
zr=zr.*(c-a+k-1.0d0).*(c-b+k-1.0d0)./(k.*(c+k-1.0d0)).*z;
zhf=zhf+zr;
end;  k=nm+1;
zhf=(1.0d0-z).^(c-a-b).*zhf;
elseif(a0 <= 1.0d0);
if(x < 0.0d0);
z1=z./(z-1.0d0);
if(c > a&b < a&b > 0.0);
a=bb;
b=aa;
end;
zc0=1.0d0./((1.0d0-z).^a);
zhf=complex(1.0d0,0.0d0);
zr0=complex(1.0d0,0.0d0);
for  k=1:500;
zr0=zr0.*(a+k-1.0d0).*(c-b+k-1.0d0)./(k.*(c+k-1.0d0)).*z1;
zhf=zhf+zr0;
if(abs(zhf-zw)< abs(zhf).*eps)break; end;
zw=zhf;
end;
zhf=zc0.*zhf;
elseif(a0 >= 0.90d0);
gm=0.0d0;
mcab=fix(c-a-b+eps.*(abs(1.0d0).*sign(c-a-b)));
if(abs(c-a-b-mcab)< eps);
m=fix(c-a-b);
[a,ga]=gamma(a,ga);
[b,gb]=gamma(b,gb);
[c,gc]=gamma(c,gc);
[dumvar1,gam]=gamma(a+m,gam);
[dumvar1,gbm]=gamma(b+m,gbm);
[a,pa]=psi(a,pa);
[b,pb]=psi(b,pb);
if(m ~= 0)gm=1.0d0; end;
for  j=1:abs(m)-1;
gm=gm.*j;
end;  j=abs(m)-1+1;
rm=1.0d0;
for  j=1:abs(m);
rm=rm.*j;
end;  j=abs(m)+1;
zf0=complex(1.0d0,0.0d0);
zr0=complex(1.0d0,0.0d0);
zr1=complex(1.0d0,0.0d0);
sp0=0.d0;
sp=0.0d0;
if(m >= 0);
zc0=gm.*gc./(gam.*gbm);
zc1=-gc.*(z-1.0d0).^m./(ga.*gb.*rm);
for  k=1:m-1;
zr0=zr0.*(a+k-1.d0).*(b+k-1.d0)./(k.*(k-m)).*(1.d0-z);
zf0=zf0+zr0;
end;  k=m-1+1;
for  k=1:m;
sp0=sp0+1.0d0./(a+k-1.0d0)+1.0./(b+k-1.0d0)-1.d0./k;
end;  k=m+1;
zf1=pa+pb+sp0+2.0d0.*el+log(1.0d0-z);
for  k=1:500;
sp=sp+(1.0d0-a)./(k.*(a+k-1.0d0))+(1.0d0-b)./(k.*(b+k-1.0d0));
sm=0.0d0;
for  j=1:m;
sm=sm+(1.0d0-a)./((j+k).*(a+j+k-1.0d0))+1.0d0./(b+j+k-1.0d0);
end;  j=m+1;
zp=pa+pb+2.0d0.*el+sp+sm+log(1.0d0-z);
zr1=zr1.*(a+m+k-1.0d0).*(b+m+k-1.0d0)./(k.*(m+k)).*(1.0d0-z);
zf1=zf1+zr1.*zp;
if(abs(zf1-zw)< abs(zf1).*eps)break; end;
zw=zf1;
end;
zhf=zf0.*zc0+zf1.*zc1;
elseif(m < 0);
m=-m;
zc0=gm.*gc./(ga.*gb.*(1.0d0-z).^m);
zc1=-(-1).^m.*gc./(gam.*gbm.*rm);
for  k=1:m-1;
zr0=zr0.*(a-m+k-1.0d0).*(b-m+k-1.0d0)./(k.*(k-m)).*(1.0d0-z);
zf0=zf0+zr0;
end;  k=m-1+1;
for  k=1:m;
sp0=sp0+1.0d0./k;
end;  k=m+1;
zf1=pa+pb-sp0+2.0d0.*el+log(1.0d0-z);
for  k=1:500;
sp=sp+(1.0d0-a)./(k.*(a+k-1.0d0))+(1.0d0-b)./(k.*(b+k-1.0d0));
sm=0.0d0;
for  j=1:m;
sm=sm+1.0d0./(j+k);
end;  j=m+1;
zp=pa+pb+2.0d0.*el+sp-sm+log(1.0d0-z);
zr1=zr1.*(a+k-1.d0).*(b+k-1.d0)./(k.*(m+k)).*(1.d0-z);
zf1=zf1+zr1.*zp;
if(abs(zf1-zw)< abs(zf1).*eps)break; end;
zw=zf1;
end;
zhf=zf0.*zc0+zf1.*zc1;
end;
else;
[a,ga]=gamma(a,ga);
[b,gb]=gamma(b,gb);
[c,gc]=gamma(c,gc);
[dumvar1,gca]=gamma(c-a,gca);
[dumvar1,gcb]=gamma(c-b,gcb);
[dumvar1,gcab]=gamma(c-a-b,gcab);
[dumvar1,gabc]=gamma(a+b-c,gabc);
zc0=gc.*gcab./(gca.*gcb);
zc1=gc.*gabc./(ga.*gb).*(1.0d0-z).^(c-a-b);
zhf=complex(0.0d0,0.0d0);
zr0=zc0;
zr1=zc1;
for  k=1:500;
zr0=zr0.*(a+k-1.d0).*(b+k-1.d0)./(k.*(a+b-c+k)).*(1.d0-z);
zr1=zr1.*(c-a+k-1.0d0).*(c-b+k-1.0d0)./(k.*(c-a-b+k)).*(1.0d0-z);
zhf=zhf+zr0+zr1;
if(abs(zhf-zw)< abs(zhf).*eps)break; end;
zw=zhf;
end;
zhf=zhf+zc0+zc1;
end;
else;
z00=complex(1.0d0,0.0d0);
if(c-a < a&c-b < b);
z00=(1.0d0-z).^(c-a-b);
a=c-a;
b=c-b;
end;
zhf=complex(1.0d0,0.d0);
zr=complex(1.0d0,0.0d0);
for  k=1:1500;
zr=zr.*(a+k-1.0d0).*(b+k-1.0d0)./(k.*(c+k-1.0d0)).*z;
zhf=zhf+zr;
if(abs(zhf-zw)<= abs(zhf).*eps)break; end;
zw=zhf;
end;
zhf=z00.*zhf;
end;
elseif(a0 > 1.0d0);
mab=fix(a-b+eps.*(abs(1.0d0).*sign(a-b)));
if(abs(a-b-mab)< eps&a0 <= 1.1d0)b=b+eps; end;
if(abs(a-b-mab)> eps);
[a,ga]=gamma(a,ga);
[b,gb]=gamma(b,gb);
[c,gc]=gamma(c,gc);
[dumvar1,gab]=gamma(a-b,gab);
[dumvar1,gba]=gamma(b-a,gba);
[dumvar1,gca]=gamma(c-a,gca);
[dumvar1,gcb]=gamma(c-b,gcb);
zc0=gc.*gba./(gca.*gb.*(-z).^a);
zc1=gc.*gab./(gcb.*ga.*(-z).^b);
zr0=zc0;
zr1=zc1;
zhf=complex(0.0d0,0.0d0);
for  k=1:500;
zr0=zr0.*(a+k-1.0d0).*(a-c+k)./((a-b+k).*k.*z);
zr1=zr1.*(b+k-1.0d0).*(b-c+k)./((b-a+k).*k.*z);
zhf=zhf+zr0+zr1;
if(abs((zhf-zw)./zhf)<= eps)break; end;
zw=zhf;
end;
zhf=zhf+zc0+zc1;
else;
if(a-b < 0.0d0);
a=bb;
b=aa;
end;
ca=c-a;
cb=c-b;
nca=fix(ca+eps.*(abs(1.0d0).*sign(ca)));
ncb=fix(cb+eps.*(abs(1.0d0).*sign(cb)));
if(abs(ca-nca)< eps|abs(cb-ncb)< eps)c=c+eps; end;
[a,ga]=gamma(a,ga);
[c,gc]=gamma(c,gc);
[dumvar1,gcb]=gamma(c-b,gcb);
[a,pa]=psi(a,pa);
[dumvar1,pca]=psi(c-a,pca);
[dumvar1,pac]=psi(a-c,pac);
mab=fix(a-b+eps);
zc0=gc./(ga.*(-z).^b);
[dumvar1,gm]=gamma(a-b,gm);
zf0=gm./gcb.*zc0;
zr=zc0;
for  k=1:mab-1;
zr=zr.*(b+k-1.0d0)./(k.*z);
t0=a-b-k;
[t0,g0]=gamma(t0,g0);
[dumvar1,gcbk]=gamma(c-b-k,gcbk);
zf0=zf0+zr.*g0./gcbk;
end;  k=mab-1+1;
if(mab == 0)zf0=complex(0.0d0,0.0d0); end;
zc1=gc./(ga.*gcb.*(-z).^a);
sp=-2.0d0.*el-pa-pca;
for  j=1:mab;
sp=sp+1.0d0./j;
end;  j=mab+1;
zp0=sp+log(-z);
sq=1.0d0;
for  j=1:mab;
sq=sq.*(b+j-1.0d0).*(b-c+j)./j;
end;  j=mab+1;
zf1=(sq.*zp0).*zc1;
zr=zc1;
rk1=1.0d0;
sj1=0.0d0;
for  k=1:10000;
zr=zr./z;
rk1=rk1.*(b+k-1.0d0).*(b-c+k)./(k.*k);
rk2=rk1;
for  j=k+1:k+mab;
rk2=rk2.*(b+j-1.0d0).*(b-c+j)./j;
end;  j=k+mab+1;
sj1=sj1+(a-1.0d0)./(k.*(a+k-1.0d0))+(a-c-1.0d0)./(k.*(a-c+k-1.0d0));
sj2=sj1;
for  j=k+1:k+mab;
sj2=sj2+1.0d0./j;
end;  j=k+mab+1;
zp=-2.0d0.*el-pa-pac+sj2-1.0d0./(k+a-c)-pi./tan(pi.*(k+a-c))+log(-z);
zf1=zf1+rk2.*zr.*zp;
ws=abs(zf1);
if(abs((ws-w0)./ws)< eps)break; end;
w0=ws;
end;
zhf=zf0+zf1;
end;
end;
a=aa;
b=bb;
if(k > 150)fprintf(1,[repmat(' ',1,1),'warning% you should check the accuracy' ' \n']); end;
%format(1x,'warning% you should check the accuracy');
return;
end
function [x,ga]=gamma(x,ga,varargin);
%     ==================================================
%     Purpose: Compute gamma function ג(x)
%     Input :  x  --- Argument of ג(x)
%(x is not equal to 0,-1,-2,תתת)
%     Output:  GA --- ג(x)
%     ==================================================
 g=zeros(1,26);
pi=3.141592653589793d0;
if(x == fix(x));
if(x > 0.0d0);
ga=1.0d0;
m1=x-1;
for  k=2:m1;
ga=ga.*k;
end;  k=m1+1;
else;
ga=1.0d+300;
end;
else;
if(abs(x)> 1.0d0);
z=abs(x);
m=fix(z);
r=1.0d0;
for  k=1:m;
r=r.*(z-k);
end;  k=m+1;
z=z-m;
else;
z=x;
end;
g(:)=[1.0d0,0.5772156649015329d0,-0.6558780715202538d0,-0.420026350340952d-1,0.1665386113822915d0,-.421977345555443d-1,-.96219715278770d-2,.72189432466630d-2,-.11651675918591d-2,-.2152416741149d-3,.1280502823882d-3,-.201348547807d-4,-.12504934821d-5,.11330272320d-5,-.2056338417d-6,.61160950d-8,.50020075d-8,-.11812746d-8,.1043427d-9,.77823d-11,-.36968d-11,.51d-12,-.206d-13,-.54d-14,.14d-14,.1d-15];
gr=g(26);
for  k=25:-1:1;
gr=gr.*z+g(k);
end;  k=1-1;
ga=1.0d0./(gr.*z);
if(abs(x)> 1.0d0);
ga=ga.*r;
if(x < 0.0d0)ga=-pi./(x.*ga.*sin(pi.*x)); end;
end;
end;
return;
end
function [x,ps]=psi(x,ps,varargin);
%     ======================================
%     Purpose: Compute Psi function
%     Input :  x  --- Argument of psi(x)
%     Output:  PS --- psi(x)
%     ======================================
xa=abs(x);
pi=3.141592653589793d0;
el=.5772156649015329d0;
s=0.0d0;
if(x == fix(x)&x <= 0.0);
ps=1.0d+300;
return;
elseif(xa == fix(xa));
n=xa;
for  k=1 :n-1;
s=s+1.0d0./k;
end;  k=n-1+1;
ps=-el+s;
elseif(xa+.5 == fix(xa+.5));
n=xa-.5;
for  k=1:n;
s=s+1.0./(2.0d0.*k-1.0d0);
end;  k=n+1;
ps=-el+2.0d0.*s-1.386294361119891d0;
else;
if(xa < 10.0);
n=10-fix(xa);
for  k=0:n-1;
s=s+1.0d0./(xa+k);
end;  k=n-1+1;
xa=xa+n;
end;
x2=1.0d0./(xa.*xa);
a1=-.8333333333333d-01;
a2=.83333333333333333d-02;
a3=-.39682539682539683d-02;
a4=.41666666666666667d-02;
a5=-.75757575757575758d-02;
a6=.21092796092796093d-01;
a7=-.83333333333333333d-01;
a8=.4432598039215686d0;
ps=log(xa)-.5d0./xa+x2.*(((((((a8.*x2+a7).*x2+a6).*x2+a5).*x2+a4).*x2+a3).*x2+a2).*x2+a1);
ps=ps-s;
end;
if(x < 0.0)ps=ps-pi.*cos(pi.*x)./sin(pi.*x)-1.0d0./x; end;
return;
end

