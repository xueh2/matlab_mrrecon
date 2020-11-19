function p = polyfit2dw(x,y,z,w,n,m)
% P= POLYFIT2D(x,y,z,w,n,m) finds the coefficients of a 
%  polynomial function of 2 variables formed from the
%  data in vectors x and y of degrees n and m, respectively,
%  that fit	the data in vector z in a weighted least-squares sense.
%  using weighting vector w
%
% The regression problem is formulated in matrix format as:
%
%	z = A*P   So that if the polynomial is cubic in  
%             x and linear in y, the problem becomes:
%
%	z = [y.*x.^3  y.*x.^2  y.*x  y  x.^3  x.^2  x  ones(length(x),1)]*
%	    [p31 p21 p11 p01 p30 p20 p10 p00]'                      
%		
%  Note that the various xy products are column vectors of length(x).
%
%  The coefficents of the output p    
%  matrix are arranged as shown:
%
%      p31 p30 
%      p21 p20 
%      p11 p10 
%      p01 p00
%
% The indices on the elements of p correspond to the 
% order of x and y associated with that element.
%
% For a solution to exist, the number of ordered 
% triples [x,y,z] must equal or exceed (n+1)*(m+1).
% Note that m or n may be zero.
%
% To evaluate the resulting polynominal function,
% use POLYVAL2D.
%
% solution minimizes weighted least squares error, e'*diag(w)*e
% where e=y-A*p and w has non-negative weights wi.

% modified 12/1/99 by Peter Kellman to add weighted least squares

% Perry W. Stout  June 29, 1995
% 4829 Rockland Way
% Fair Oaks, CA  95628
% (916) 966-0236
% Based on the Matlab function polyfit.

global a_vandermonde; % make this variable global to reduce repeated calculation

if any((size(x) ~= size(y)) | (size(z) ~= size(y)))
	error('X, Y,and Z vectors must be the same size')
end

x = x(:); y = y(:); z= z(:);  w=w(:);% Switches vectors to columns--matrices, too
w=double(w);

if length(x) < (n+1)*(m+1)
 error('Number of points must equal or exceed order of polynomial function.')
end

n = n + 1;
m = m + 1; % Increments n and m to equal row, col numbers of p.

if (exist('a_vandermonde')==1 & min(size(a_vandermonde)==[max(size(x)),n*m]));
else
   a_vandermonde = zeros(max(size(x)),n*m);
	% Construct the extended Vandermonde matrix, containing all xy products
	for i1= 1:m
  		for j1=1:n
	   	a_vandermonde(:,j1+(i1-1)*n) = (x.^(n-j1)).*(y.^(m-i1));
   	end
   end
end

% weighting added by Kellman 12/30/99
if 0
    sqrtW=diag(sqrt(w));
	aa=sqrtW*a_vandermonde;
	zz=sqrtW*z;
end
sqrtW=sqrt(w);
aa=zeros(size(a_vandermonde));
for i=1:size(a_vandermonde,2)
   aa(:,i)=sqrtW.*a_vandermonde(:,i);
end
zz=sqrtW.*z;


% method used in version by Stout based on MATLAB \ for least squares solution
%p1 = (aa\zz); % p1=(a\z) is unweighted least squares solution
% Reform p as a matrix.

    % modified by Kellman 12/30/99 to use factorization as done
    % in MATLAB function polyfit
    [Q R]=qr(aa,0);
    p1=R\(Q'*zz);

% p1=pinv(aa)*zz; % using pseudo-inverse function, pinv()

p=[];
for i1=1:m
	p=[p, p1((n*(i1-1)+1):(n*i1))];
end

