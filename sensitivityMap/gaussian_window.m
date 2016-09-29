function w=gaussian_window(n, sigma, sflag, circularsupport);
% function w=gaussian_window(n, sigma,'sflag', circularsupport);
% 
%   calculated 1-d or 2-d Gaussian window truncated at sigma
%
%   [n] for 1-d window
%   [n1 n2] for 2-d window (n1 rows, n2 columns)
%   sigma (default = 1) is the (truncation) value of window at edges
%   sflag is either 'symmetric' (default) or 'periodic'
%   circularsupport (2-d case) is a {0,1} flag for region of support mask (default=0, rectangular support)
%
%   default usage:
%
%   gaussian_window(N) returns the N-point symmetric gaussian window in a column vector with sigma=1
%   gaussian_window([N1 N2]) returns a (symmetric) 2-d window with N1 rows and N2 columns with sigma=1
%
%   1-d Gaussian window characteristics (same for 2-d w/ rectangular suport):
%   sigma    FWHM    max sidelobe amplitude
%     1.0     1.3          0.146
%     1.5     1.4          0.076
%     2.0     1.7          0.024
%     2.5     1.95         0.0063
%   uniform window
%     n/a     1.21         0.21

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  NIH NHLBI                          *
%     ***************************************

% set defaults for sigma and sflag
if nargin==1; sigma=1; sflag='symmetric'; circularsupport=0; end
if nargin==2; sflag='symmetric'; circularsupport=0; end
if nargin==3; circularsupport=0; end

a=sigma^2/2; % note sigma is in the numerator since x&y are scaled between [-1,1]
if length(n)==1 % 1-d window case
	switch sflag
	case 'symmetric'
        x=linspace(-1,1,n)';
        w=exp(-a*x.^2);
	case 'periodic'
        x=linspace(-1,1,n+1)';
        w=exp(-a*x.^2);
        w=w(1:end-1);
	end
elseif length(n)==2; % 2-d window case
	switch sflag
	case 'symmetric'
		x=linspace(-1,1,n(2));
		y=linspace(-1,1,n(1));
		[X,Y]=meshgrid(x,y);
		rsquared=((X).^2+(Y).^2);
        w=exp(-a*rsquared);
	case 'periodic'
		x=linspace(-1,1,n(2)+1);
		y=linspace(-1,1,n(1)+1);
		[X,Y]=meshgrid(x,y);
		rsquared=((X).^2+(Y).^2);
        w=exp(-a*rsquared);
        w=w(1:end-1,1:end-1);
	end
    if circularsupport==1;
        e=elliptical_support_region(n(1),n(2));
        w=w.*e;
    end
else
    disp('Only 1-d and 2-d windows calculated using gaussian_window')
    w=[];
end


function e=elliptical_support_region(rows,cols);
% function e=elliptical_support_region(ysize,xsize);
%
% function to compute mask for elliptical support region
% circumscribed by rectangle with dimensions rows x cols


if nargin==1; cols=rows; end
x=linspace(-(cols-1)/cols,(cols-1)/cols,cols);
y=linspace(-(rows-1)/rows,(rows-1)/rows,rows);
[X,Y]=meshgrid(x,y);
r=sqrt(X.^2+Y.^2);
e=ones(rows,cols);
e(r>1)=0;

return

