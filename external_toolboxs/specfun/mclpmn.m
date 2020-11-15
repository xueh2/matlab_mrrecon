function mclpmn
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ============================================================
%       Purpose: This program computes the associated Legendre
%                functions Pmn(z)and their derivatives Pmn'(z)for
%                a complex argument using subroutine CLPMN
%       Input :  x --- Real part of z
%                y --- Imaginary part of z
%                m --- Order of Pmn(z),  m = 0,1,2,...,n
%                n --- Degree of Pmn(z), n = 0,1,2,...,N
%       Output:  CPM(m,n)--- Pmn(z)
%                CPD(m,n)--- Pmn'(z)
%       Examples:
%                n = 5, x = 0.5, y = 0.2
%       m     Re[Pmn(z)]Im[Pmn(z)]Re[Pmn'(z)]Im[Pmn'(z)]
%      -------------------------------------------------------------
%       0    .252594D+00  -.530293D+00   -.347606D+01  -.194250D+01
%       1    .333071D+01   .135206D+01    .117643D+02  -.144329D+02
%       2   -.102769D+02   .125759D+02    .765713D+02   .598500D+02
%       3   -.572879D+02  -.522744D+02   -.343414D+03   .147389D+03
%       4    .335711D+03  -.389151D+02   -.226328D+03  -.737100D+03
%       5   -.461125D+03   .329122D+03    .187180D+04   .160494D+02
%                n = 5, x = 2.5, y = 1.0
%       m     Re[Pmn(z)]Im[Pmn(z)]Re[Pmn'(z)]Im[Pmn'(z)]
%      -------------------------------------------------------------
%       0   -.429395D+03   .900336D+03   -.350391D+02   .193594D+04
%       1   -.216303D+04   .446358D+04   -.208935D+03   .964685D+04
%       2   -.883477D+04   .174005D+05   -.123703D+04   .381938D+05
%       3   -.273211D+05   .499684D+05   -.568080D+04   .112614D+06
%       4   -.565523D+05   .938503D+05   -.167147D+05   .219713D+06
%       5   -.584268D+05   .863328D+05   -.233002D+05   .212595D+06
%       ============================================================
m=[];n=[];x=[];y=[];cpm=[];cpd=[];
 cpm=zeros(40+1,40+1);
cpd=zeros(40+1,40+1);
fprintf(1,'%s \n','  please enter m, n, x and y ');
%        READ(*,*)M,N,X,Y
m=0;
n=5;
x=.5;
y=.2;
fprintf(1,[repmat(' ',1,1),'m =','%2g',', ','n =','%2g',', ','x =','%5.1g',', ','y =','%5.1g' ' \n'],m,n,x,y);
[dumvar1,m,n,x,y,cpm,cpd]=clpmn(40,m,n,x,y,cpm,cpd);
fprintf(1,'%s ','   m   n    re[pmn(z)]im[pmn(z)]');fprintf(1,'%s \n','re[pmn''(z)]im[pmn''(z)]');
fprintf(1,'%s ',' -----------------------------------');fprintf(1,'%s \n','-------------------------------');
for  j=0:n;
fprintf(1,[repmat(' ',1,1),repmat('%4g',1,2),repmat(' ',1,1),repmat('%14.6g',1,2),repmat(' ',1,1),repmat('%14.6g',1,2) ' \n'],m,j,cpm(m+1,j+1),cpd(m+1,j+1));
end;  j=n+1;
%format(1x,2i4,1x,2d14.6,1x,2d14.6);
%format(1x,',i2,', ',',i2,', ',',f5.1, ', ',',f5.1);
end
function [mm,m,n,x,y,cpm,cpd]=clpmn(mm,m,n,x,y,cpm,cpd,varargin);
%       =========================================================
%       Purpose: Compute the associated Legendre functions Pmn(z)
%                and their derivatives Pmn'(z)for a complex
%                argument
%       Input :  x  --- Real part of z
%                y  --- Imaginary part of z
%                m  --- Order of Pmn(z),  m = 0,1,2,...,n
%                n  --- Degree of Pmn(z), n = 0,1,2,...,N
%                mm --- Physical dimension of CPM and CPD
%       Output:  CPM(m,n)--- Pmn(z)
%                CPD(m,n)--- Pmn'(z)
%       =========================================================
z=complex(x,y);
for  i=0:n;
for  j=0:m;
cpm(j+1,i+1)=complex(0.0d0,0.0d0);
cpd(j+1,i+1)=complex(0.0d0,0.0d0);
end;  j=fix(m)+1;
end;  i=fix(n)+1;
cpm(0+1,0+1)=complex(1.0d0,0.0d0);
if(abs(x)== 1.0d0&y == 0.0d0);
for  i=1:n;
cpm(0+1,i+1)=x.^i;
cpd(0+1,i+1)=0.5d0.*i.*(i).*x.^(i);
end;  i=fix(n)+1;
for  j=1:n;
for  i=1:m;
if(i == 1);
cpd(i+1,j+1)=complex(1.0d+300,0.0d0);
elseif(i == 2);
cpd(i+1,j+1)=-0.25d0.*(j+2).*(j).*j.*(j-1).*x.^(j);
end;
end;  i=fix(m)+1;
end;  j=fix(n)+1;
return;
end;
ls=1;
if(abs(z)> 1.0d0)ls=-1; end;
zq=sqrt(ls.*(1.0d0-z.*z));
zs=ls.*(1.0d0-z.*z);
for  i=1:m;
cpm(i+1,i+1)=-ls.*(2.0d0.*i-1.0d0).*zq.*cpm(i-1+1,i-1+1);
end;  i=fix(m)+1;
for  i=0:m;
cpm(i+1,i+1+1)=(2.0d0.*i+1.0d0).*z.*cpm(i+1,i+1);
end;  i=fix(m)+1;
for  i=0:m;
for  j=i+2:n;
cpm(i+1,j+1)=((2.0d0.*j-1.0d0).*z.*cpm(i+1,j-1+1)-(i+j- 1.0d0).*cpm(i+1,j-2+1))./(j-i);
end;  j=fix(n)+1;
end;  i=fix(m)+1;
cpd(0+1,0+1)=complex(0.0d0,0.0d0);
for  j=1:n;
cpd(0+1,j+1)=ls.*j.*(cpm(0+1,j-1+1)-z.*cpm(0+1,j+1))./zs;
end;  j=fix(n)+1;
for  i=1:m;
for  j=i:n;
cpd(i+1,j+1)=ls.*i.*z.*cpm(i+1,j+1)./zs+(j+i).*(j-i+1.0d0)./zq.*cpm(i-1+1,j+1);
end;  j=fix(n)+1;
end;  i=fix(m)+1;
return;
end

