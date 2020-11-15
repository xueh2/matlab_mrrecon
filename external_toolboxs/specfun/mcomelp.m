function mcomelp
%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%

%       ===================================================
%       Purpose: This program computes complete elliptic
%                integrals K(k)and E(k)using subroutine
%                COMELP
%       Input  : K  --- Modulus k(0 ó k ó 1)
%       Output : CK --- K(k)
%                CE --- E(k)
%       Example:
%                  k         K(k)E(K)
%                ---------------------------------
%                 .00      1.570796      1.570796
%                 .25      1.596242      1.545957
%                 .50      1.685750      1.467462
%                 .75      1.910990      1.318472
%                1.00       ì            1.000000
%       ===================================================
hk=[];ck=[];ce=[];
  hk=0;
 ck=0;
 ce=0;
fprintf(1,'%s \n','please enter the modulus k ');
%        READ(*,*)HK
hk=.5;
fprintf(1,'%s \n','    k         k(k)e(k)');
fprintf(1,'%s \n','  ---------------------------------');
[hk,ck,ce]=comelp(hk,ck,ce);
if(hk ~= 1.0)fprintf(1,[repmat(' ',1,2),'%5.2g',repmat('%14.6g',1,2) ' \n'],hk,ck,ce); end;
if(hk == 1.0)fprintf(1,[repmat(' ',1,2),'%5.2g',repmat(' ',1,7),'ì',repmat(' ',1,6),'%14.6g' ' \n'],hk,ce); end;
%format(2x,f5.2,2f14.6);
%format(2x,f5.2,7x,'ì',6x,f14.6);
end
function [hk,ck,ce]=comelp(hk,ck,ce,varargin);
%       ==================================================
%       Purpose: Compute complete elliptic integrals K(k)
%                and E(k)
%       Input  : K  --- Modulus k(0 ó k ó 1)
%       Output : CK --- K(k)
%                CE --- E(k)
%       ==================================================
pk=1.0d0-hk.*hk;
if(hk == 1.0);
ck=1.0d+300;
ce=1.0d0;
else;
ak=(((.01451196212d0.*pk+.03742563713d0).*pk+.03590092383d0).*pk+.09666344259d0).*pk+ 1.38629436112d0;
bk=(((.00441787012d0.*pk+.03328355346d0).*pk+.06880248576d0).*pk+.12498593597d0).*pk+.5d0;
ck=ak-bk.*log(pk);
ae=(((.01736506451d0.*pk+.04757383546d0).*pk+.0626060122d0).*pk+.44325141463d0).*pk+1.0d0;
be=(((.00526449639d0.*pk+.04069697526d0).*pk+.09200180037d0).*pk+.2499836831d0).*pk;
ce=ae-be.*log(pk);
end;
return;
end

