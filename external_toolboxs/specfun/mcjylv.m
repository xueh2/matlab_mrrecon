function mcjylv
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ========================================================
%       Purpose: This program computes Bessel functions Jv(z)
%                and Yv(z)and their derivatives with a large
%                order and complex argument using subroutine
%                CJYLV
%       Input:   v --- Order of Jv(z)and Yv(z)
%                z --- Complex argument
%       Output:  CBJV --- Jv(z)
%                CDJV --- Jv'(z)
%                CBYV --- Yv(z)
%                CDYV --- Yv'(z)
%       Examples:
%                v = 100.00,    z = 4.00 + 2.00 i
%                Jv(z)= -.6444792518-123 + .6619157435-123 i
%                Jv'(z)= -.6251103777-122 + .1967638668-121 i
%                Yv(z)=  .2403065353+121 + .2472039414+121 i
%                Yv'(z)= -.7275814786+122 - .2533588851+122 i
%                v =100.5,     z = 4.00 + 2.00 i
%                Jv(z)= -.1161315754-123 + .7390127781-124 i
%                Jv'(z)= -.1588519437-122 + .2652227059-122 i
%                Yv(z)=  .1941381412+122 + .1237578195+122 i
%                Yv'(z)= -.5143285247+123 - .5320026773+122 i
%       ========================================================
v=[];z=[];cbjv=[];cdjv=[];cbyv=[];cdyv=[];
fprintf(1,'%s \n','please enter v,x and y(z = x+iy)');
%        READ(*,*)V,X,Y
v=100.0;
x=4.0;
y=2.0;
fprintf(1,[repmat(' ',1,8),'v = ','%6.2g',',    ','z =','%7.2g',' + i ','%7.2g' ' \n'],v,x,y);
z=complex(x,y);
[v,z,cbjv,cdjv,cbyv,cdyv]=cjylv(v,z,cbjv,cdjv,cbyv,cdyv);
fprintf(1,'%0.15g \n');
fprintf(1,[repmat(' ',1,8),'jv(z)=','%17.10g',' + i','%17.10g' ' \n'],cbjv);
fprintf(1,[repmat(' ',1,8),'jv''(z)=','%17.10g',' + i','%17.10g' ' \n'],cdjv);
fprintf(1,'%0.15g \n');
fprintf(1,[repmat(' ',1,8),'yv(z)=','%17.10g',' + i','%17.10g' ' \n'],cbyv);
fprintf(1,[repmat(' ',1,8),'yv''(z)=','%17.10g',' + i','%17.10g' ' \n'],cdyv);
%format(8x,'v = ',f6.2,',    ',',f7.2,' + i ',f7.2);
%format(8x,',d17.10,' + i',d17.10);
%format(8x,',d17.10,' + i',d17.10);
%format(8x,',d17.10,' + i',d17.10);
%format(8x,',d17.10,' + i',d17.10);
end
function [v,z,cbjv,cdjv,cbyv,cdyv]=cjylv(v,z,cbjv,cdjv,cbyv,cdyv,varargin);
%       ===================================================
%       Purpose: Compute Bessel functions Jv(z)and Yv(z)
%                and their derivatives with a complex
%                argument and a large order
%       Input:   v --- Order of Jv(z)and Yv(z)
%                z --- Complex argument
%       Output:  CBJV --- Jv(z)
%                CDJV --- Jv'(z)
%                CBYV --- Yv(z)
%                CDYV --- Yv'(z)
%       Routine called:
%                CJK to compute the expansion coefficients
%       ===================================================
km=[];a=[];
 cf=zeros(1,12);
a=zeros(1,91);
km=12;
[km,a]=cjk(km,a);
pi=3.141592653589793d0;
for  l=1:-1:0;
v0=v-l;
cws=sqrt(1.0d0-(z./v0).*(z./v0));
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
csj=complex(1.0d0,0.0d0);
for  k=1:km;
csj=csj+cf(k).*vr.^k;
end;  k=km+1;
cbjv=sqrt(ct./(2.0d0.*pi.*v0)).*exp(v0.*ceta).*csj;
if(l == 1)cfj=cbjv; end;
csy=complex(1.0d0,0.0d0);
for  k=1:km;
csy=csy+(-1).^k.*cf(k).*vr.^k;
end;  k=km+1;
cbyv=-sqrt(2.0d0.*ct./(pi.*v0)).*exp(-v0.*ceta).*csy;
if(l == 1)cfy=cbyv; end;
end;  l=0-1;
cdjv=-v./z.*cbjv+cfj;
cdyv=-v./z.*cbyv+cfy;
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

