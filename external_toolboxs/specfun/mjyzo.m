function mjyzo
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ==========================================================
%       Purpose: This program computes the zeros of Bessel
%                functions Jn(x), Yn(x), and their derivatives
%                using subroutine JYZO
%       Input :  n --- Order of Bessel functions(n ó 100)
%                NT --- Number of zeros
%       Output:  RJ0(m)--- m-th zero of Jn(x),  m=1,2,...,NT
%                RJ1(m)--- m-th zero of Jn'(x), m=1,2,...,NT
%                RY0(m)--- m-th zero of Yn(x),  m=1,2,...,NT
%                RY1(m)--- m-th zero of Yn'(x), m=1,2,...,NT
%       Example: n = 1, NT =5
%      Zeros of Bessel funcions Jn(x), Yn(x)and their derivatives
%(n = 1)
%       m       jnm           j'nm          ynm           y'nm
%      -----------------------------------------------------------
%       1     3.8317060     1.8411838     2.1971413     3.6830229
%       2     7.0155867     5.3314428     5.4296810     6.9415000
%       3    10.1734681     8.5363164     8.5960059    10.1234047
%       4    13.3236919    11.7060049    11.7491548    13.2857582
%       5    16.4706301    14.8635886    14.8974421    16.4400580
%       ==========================================================
n=[];nt=[];rj0=[];rj1=[];ry0=[];ry1=[];
 rj0=zeros(1,101);
rj1=zeros(1,101);
ry0=zeros(1,101);
ry1=zeros(1,101);
fprintf(1,'%s \n','please enter n and nt ');
%        READ(*,*)N,NT
n=1;
nt=5;
fprintf(1,'%0.15g \n');
[n,nt,rj0,rj1,ry0,ry1]=jyzo(n,nt,rj0,rj1,ry0,ry1);
fprintf(1,[repmat(' ',1,2),'zeros of bessel funcions jn(x), yn(x)',' and their derivatives' ' \n']);
fprintf(1,[repmat(' ',1,30),'(n =','%2g',')' ' \n'],n);
fprintf(1,'%s ','  m       jnm           j''nm          ynm');fprintf(1,'%s \n', '           y''nm');
fprintf(1,'%s ',' ----------------------------------------');fprintf(1,'%s \n', '-------------------');
for  m=1:nt;
fprintf(1,[repmat(' ',1,1),'%3g',repmat('%14.7g',1,4) ' \n'],m,rj0(m),rj1(m),ry0(m),ry1(m));
end;  m=nt+1;
%format(2x,'zeros of bessel funcions jn(x), yn(x)',' and their derivatives');
%format(30x,',i2,')');
%format(1x,i3,4f14.7);
end
function [n,nt,rj0,rj1,ry0,ry1]=jyzo(n,nt,rj0,rj1,ry0,ry1,varargin);
%       ======================================================
%       Purpose: Compute the zeros of Bessel functions Jn(x),
%                Yn(x), and their derivatives
%       Input :  n  --- Order of Bessel functions(n ó 101)
%                NT --- Number of zeros(roots)
%       Output:  RJ0(L)--- L-th zero of Jn(x),  L=1,2,...,NT
%                RJ1(L)--- L-th zero of Jn'(x), L=1,2,...,NT
%                RY0(L)--- L-th zero of Yn(x),  L=1,2,...,NT
%                RY1(L)--- L-th zero of Yn'(x), L=1,2,...,NT
%       Routine called: JYNDD for computing Jn(x), Yn(x), and
%                       their first and second derivatives
%       ======================================================
x=[];bjn=[];djn=[];fjn=[];byn=[];dyn=[];fyn=[];
if(n <= 20);
x=2.82141+1.15859.*fix(n);
else;
x=fix(n)+1.85576.*fix(n).^0.33333+1.03315./fix(n).^0.33333;
end;
l=0;
while (1);
while (1);
x0=x;
[n,x,bjn,djn,fjn,byn,dyn,fyn]=jyndd(fix(n),x,bjn,djn,fjn,byn,dyn,fyn);
x=x-bjn./djn;
if(~(abs(x-x0)> 1.0d-9))break; end;
end;
l=l+1;
rj0(l)=x;
x=x+3.1416+(0.0972+0.0679.*fix(n)-0.000354.*fix(n).^2)./l;
if(~(l < nt))break; end;
end;
if(n <= 20);
x=0.961587+1.07703.*fix(n);
else;
x=fix(n)+0.80861.*fix(n).^0.33333+0.07249./fix(n).^0.33333;
end;
if(n == 0)x=3.8317; end;
l=0;
while (1);
while (1);
x0=x;
[n,x,bjn,djn,fjn,byn,dyn,fyn]=jyndd(fix(n),x,bjn,djn,fjn,byn,dyn,fyn);
x=x-djn./fjn;
if(~(abs(x-x0)> 1.0d-9))break; end;
end;
l=l+1;
rj1(l)=x;
x=x+3.1416+(0.4955+0.0915.*fix(n)-0.000435.*fix(n).^2)./l;
if(~(l < nt))break; end;
end;
if(n <= 20);
x=1.19477+1.08933.*fix(n);
else;
x=fix(n)+0.93158.*fix(n).^0.33333+0.26035./fix(n).^0.33333;
end;
l=0;
while (1);
while (1);
x0=x;
[n,x,bjn,djn,fjn,byn,dyn,fyn]=jyndd(fix(n),x,bjn,djn,fjn,byn,dyn,fyn);
x=x-byn./dyn;
if(~(abs(x-x0)> 1.0d-9))break; end;
end;
l=l+1;
ry0(l)=x;
x=x+3.1416+(0.312+0.0852.*fix(n)-0.000403.*fix(n).^2)./l;
if(~(l < nt))break; end;
end;
if(n <= 20);
x=2.67257+1.16099.*fix(n);
else;
x=fix(n)+1.8211.*fix(n).^0.33333+0.94001./fix(n).^0.33333;
end;
l=0;
while (1);
while (1);
x0=x;
[n,x,bjn,djn,fjn,byn,dyn,fyn]=jyndd(fix(n),x,bjn,djn,fjn,byn,dyn,fyn);
x=x-dyn./fyn;
if(~(abs(x-x0)> 1.0d-9))break; end;
end;
l=l+1;
ry1(l)=x;
x=x+3.1416+(0.197+0.0643.*fix(n)-0.000286.*fix(n).^2)./l;
if(~(l < nt))break; end;
end;
return;
end
function [n,x,bjn,djn,fjn,byn,dyn,fyn]=jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn,varargin);
%       ===========================================================
%       Purpose: Compute Bessel functions Jn(x)and Yn(x), and
%                their first and second derivatives
%       Input:   x   ---  Argument of Jn(x)and Yn(x,x > 0)
%                n   ---  Order of Jn(x)and Yn(x)
%       Output:  BJN ---  Jn(x)
%                DJN ---  Jn'(x)
%                FJN ---  Jn"(x)
%                BYN ---  Yn(x)
%                DYN ---  Yn'(x)
%                FYN ---  Yn"(x)
%       ===========================================================
 bj=zeros(1,102);
by=zeros(1,102);
for  nt=1:900;
mt=fix(0.5.*log10(6.28.*nt)-nt.*log10(1.36.*abs(x)./nt));
if(mt > 20)break; end;
end;
m=nt;
bs=0.0d0;
f0=0.0d0;
f1=1.0d-35;
su=0.0d0;
for  k=m:-1:0;
f=2.0d0.*(k+1.0d0).*f1./x-f0;
if(k <= n+1)bj(k+1)=f; end;
if(k == 2.*fix(k./2));
bs=bs+2.0d0.*f;
if(k ~= 0)su=su+(-1).^(k./2).*f./k; end;
end;
f0=f1;
f1=f;
end;  k=0-1;
for  k=0:n+1;
bj(k+1)=bj(k+1)./(bs-f);
end;  k=fix(n)+1+1;
bjn=bj(fix(n)+1);
ec=0.5772156649015329d0;
e0=0.3183098861837907d0;
s1=2.0d0.*e0.*(log(x./2.0d0)+ec).*bj(1);
f0=s1-8.0d0.*e0.*su./(bs-f);
f1=(bj(2).*f0-2.0d0.*e0./x)./bj(1);
by(1)=f0;
by(2)=f1;
for  k=2:n+1;
f=2.0d0.*(k-1.0d0).*f1./x-f0;
by(k+1)=f;
f0=f1;
f1=f;
end;  k=fix(n)+1+1;
byn=by(fix(n)+1);
djn=-bj(fix(n)+2)+fix(n).*bj(fix(n)+1)./x;
dyn=-by(fix(n)+2)+fix(n).*by(fix(n)+1)./x;
fjn=(fix(n).*fix(n)./(x.*x)-1.0d0).*bjn-djn./x;
fyn=(fix(n).*fix(n)./(x.*x)-1.0d0).*byn-dyn./x;
return;
end

