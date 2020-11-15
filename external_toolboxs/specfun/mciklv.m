function mciklv
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       =========================================================
%       Purpose: This program computes modified Bessel functions
%                Iv(z)and Kv(z)and their derivatives for a
%                large order and a complex argument using
%                subroutine CIKLV
%       Input:   v --- Order of Iv(z)and Kv(z)
%                z --- Complex argument
%       Output:  CBIV --- Iv(z)
%                CDIV --- Iv'(z)
%                CBKV --- Kv(z)
%                CDKV --- Kv'(z)
%       Examples:
%                v =100.00,    z =   4.00 + i   2.00
%       Iv(z)= -.7373606617-123 + .6461109082-123 i
%       Iv'(z)= -.8307094243-122 + .2030132500-121 i
%       Kv(z)= -.3836166007+121 - .3356017795+121 i
%       Kv'(z)=  .1103271276+123 + .2886519240+122 i
%                v =100.50,    z =   4.00 + i   2.00
%       Iv(z)= -.1289940051-123 + .6845756182-124 i
%       Iv'(z)= -.1907996261-122 + .2672465997-122 i
%       Kv(z)= -.3008779281+122 - .1593719779+122 i
%       Kv'(z)=  .7653781978+123 + .1857772148+122 i
%       =========================================================
v=[];z=[];cbiv=[];cdiv=[];cbkv=[];cdkv=[];
fprintf(1,'%s \n','please enter v,x,y(z = x+iy)');
%        READ(*,*)V,X,Y
v=100.00;
x=4.0;
y=2.0;
fprintf(1,[repmat(' ',1,8),'v =','%6.2g',',    ','z =','%7.2g',' + i','%7.2g' ' \n'],v,x,y);
z=complex(x,y);
[v,z,cbiv,cdiv,cbkv,cdkv]=ciklv(v,z,cbiv,cdiv,cbkv,cdkv);
fprintf(1,'%0.15g \n');
fprintf(1,[repmat(' ',1,8),'iv(z)=','%17.10g',' + i ','%17.10g' ' \n'],cbiv);
fprintf(1,[repmat(' ',1,8),'iv''(z)=','%17.10g',' + i ','%17.10g' ' \n'],cdiv);
fprintf(1,'%0.15g \n');
fprintf(1,[repmat(' ',1,8),'kv(z)=','%17.10g',' + i ','%17.10g' ' \n'],cbkv);
fprintf(1,[repmat(' ',1,8),'kv''(z)=','%17.10g',' + i ','%17.10g' ' \n'],cdkv);
%format(8x,',f6.2,',    ',',f7.2,' + i',f7.2);
%format(8x,',d17.10,' + i ',d17.10);
%format(8x,',d17.10,' + i ',d17.10);
%format(8x,',d17.10,' + i ',d17.10);
%format(8x,',d17.10,' + i ',d17.10);
end
function [v,z,cbiv,cdiv,cbkv,cdkv]=ciklv(v,z,cbiv,cdiv,cbkv,cdkv,varargin);
%       =====================================================
%       Purpose: Compute modified Bessel functions Iv(z)and
%                Kv(z)and their derivatives with a complex
%                argument and a large order
%       Input:   v --- Order of Iv(z)and Kv(z)
%                z --- Complex argument
%       Output:  CBIV --- Iv(z)
%                CDIV --- Iv'(z)
%                CBKV --- Kv(z)
%                CDKV --- Kv'(z)
%       Routine called:
%                CJK to compute the expansion coefficients
%       ====================================================
km=[];a=[];
 cf=zeros(1,12);
a=zeros(1,91);
pi=3.141592653589793d0;
km=12;
[km,a]=cjk(km,a);
for  l=1:-1:0;
v0=v-l;
cws=sqrt(1.0d0+(z./v0).*(z./v0));
ceta=cws+log(z./v0./(1.0d0+cws));
ct=1.0d0./cws;
ct2=ct.*ct;
for  k=1:km;
l0=k.*(k+1)./2+1;
lf=l0+k;
cf(k)=a(lf);
for  i=lf-1:-1:l0;
cf(k)=cf(k).*ct2+a(i);
end;  i=l0-1;
cf(k)=cf(k).*ct.^k;
end;  k=km+1;
vr=1.0d0./v0;
csi=complex(1.0d0,0.0d0);
for  k=1:km;
csi=csi+cf(k).*vr.^k;
end;  k=km+1;
cbiv=sqrt(ct./(2.0d0.*pi.*v0)).*exp(v0.*ceta).*csi;
if(l == 1)cfi=cbiv; end;
csk=complex(1.0d0,0.0d0);
for  k=1:km;
csk=csk+(-1).^k.*cf(k).*vr.^k;
end;  k=km+1;
cbkv=sqrt(pi.*ct./(2.0d0.*v0)).*exp(-v0.*ceta).*csk;
if(l == 1)cfk=cbkv; end;
end;  l=0-1;
cdiv=cfi-v./z.*cbiv;
cdkv=-cfk-v./z.*cbkv;
return;
end
function [km,a]=cjk(km,a,varargin);
%       ========================================================
%       Purpose: Compute the expansion coefficients for the
%                asymptotic expansion of Bessel functions
%                with large orders
%       Input :  Km   --- Maximum k
%       Output:  A(L)--- Cj(k)where j and k are related to L
%                         by L=j+1+[k*(k+1)]/2; j,k=0,1,...,Km
%       ========================================================
a(1)=1.0d0;
f0=1.0d0;
g0=1.0d0;
for  k=0:km-1;
l1=(k+1).*(k+2)./2+1;
l2=(k+1).*(k+2)./2+k+2;
f=(0.5d0.*k+0.125d0./(k+1)).*f0;
g=-(1.5d0.*k+0.625d0./(3.0.*(k+1.0d0))).*g0;
a(l1)=f;
a(l2)=g;
f0=f;
g0=g;
end;  k=fix(km)-1+1;
for  k=1:km-1;
for  j=1:k;
l3=k.*(k+1)./2+j+1;
l4=(k+1).*(k+2)./2+j+1;
a(l4)=(j+0.5d0.*k+0.125d0./(2.0.*j+k+1.0)).*a(l3)-(j+0.5d0.*k-1.0+0.625d0./(2.0.*j+k+1.0)).*a(l3-1);
end;  j=k+1;
end;  k=fix(km)-1+1;
return;
end

