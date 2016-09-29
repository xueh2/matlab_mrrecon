function [x, prevXs] = pocsSPIRiT(data, GOP, nIter, x0, wavWeight, show, stopThres, centre, width, continuationStep, wavThresRatio)
%
%
% res = pocsSPIRIT(y, GOP, nIter, x0, wavWeight, show)
%
% Implementation of the Cartesian, POCS l1-SPIRiT reconstruction
%
% Input:
%		y - Undersampled k-space data. Make sure that empty entries are zero
%			or else they will not be filled.
%		GOP -	the SPIRiT kernel operator obtained by calibration
%		nIter -	Maximum number of iterations
%		x0 - initial guess
%		wavWeight - wavlet threshodling parameter
% 		show - >1 to display progress (slower)
%
% Outputs:
%		res - Full k-space matrix
%
% (c) Michael Lustig 2007
%
if nargin<7
    stopThres = 0.0022;
end

if nargin<9
    centre = 250;
    width = 500;
end

if nargin<11
    continuationStep = 10;
    wavThresRatio = 1.5;
end

if nargout>1
    prevXs = cell(nIter, 1);
end

% if no l1 penalt then skip wavelet thresholding.
if wavWeight==0

	mask = (data==0);

	x = x0;
    prevX = x;
    if nargout>1 prevXs{1} = prevX; end
    
	for n=1:nIter
        
        prevX = x;
		if nargout>1 prevXs{n} = prevX; end
        
        tmpx =(x + GOP*x).*(mask); % Apply (G-I)x + x
		x = tmpx + data; % fix the data        
		if show
            if ( n == 1 ) 
                h=-1 
            end
            h = plotKSpaceInFigure(x, [1 1 1], centre, width, 1, 1/24, h);
            norm(x(:)-prevX(:))./norm(x(:))
        end
        err = norm(x(:)-prevX(:))./norm(x(:));
        if ( err < stopThres )
            disp(['Pocs SPIRiT stops at step' num2str(n)]);
            break;
        end
	end

else
    
    % find the closest diadic size for the images
	[sx,sy,nc] = size(data);
	ssx = 2^ceil(log2(sx)); 
	ssy = 2^ceil(log2(sy));
	ss = max(ssx, ssy);
    
	% W = Wavelet('Daubechies',4,4);
	% W = Wavelet('Haar',2,3);
    
    W = UndecimatedWavelet(2, 'db1', 'ppd');

	mask = (data==0);
	x = x0;
	prevX = x0;
    if nargout>1 prevXs {1} = prevX; end
    
    % Ding's suggestion
    alpha = 0.5;
    
	for n=1:nIter
        
        tic
        
        prevX = x;
        if nargout>1 prevXs{n} = prevX; end
        
		x = (x + GOP*x ).*(mask) + data; % Apply (G-I)*x + x        
        x = alpha*prevX + (1-alpha)*x;
        
        % apply wavelet thresholding
        X = ifft2c(x); % goto image domain
 		%X= zpad(X,ss,ss,nc); % zpad to the closest diadic 
		X = W*(X); % apply wavelet
		
        % X = softThresh(X,wavWeight); % threshold ( joint sparsity)
        X = W.softThresh(X,wavWeight); % threshold ( joint sparsity)
		
        X = W'*(X); % get back the image
 		% X = crop(X,sx,sy,nc); % return to the original size
		xx = fft2c(X); % go back to k-space
		x = xx.*mask + data; % fix the data
		
        timeUsed = toc;
		if show
            if ( n == 1 ) 
                h=-1 
            end
            h = plotKSpaceInFigure(x, [1 1 1], centre, width, 1, 1/24, h);
            disp([ 'Time used for iter ' num2str(n) ' : ' num2str(timeUsed) 's : delta - ' num2str(norm(x(:)-prevX(:))./norm(x(:))) ])
        end
        
        err = norm(x(:)-prevX(:))./norm(x(:));
        if ( err < stopThres )
            disp(['Pocs SPIRiT stops at step' num2str(n)]);
            break;
        end
        
        % add the continuation
        if ( mod(n, continuationStep) == 0 )
            wavWeight = wavWeight / wavThresRatio;
        end
    end
end

if nargout>1 
    prevXs = prevXs(1:n); 
end

function x = softThresh(y,t)
% apply joint sparsity soft-thresholding 
absy = sqrt(sum(abs(y).^2,3));
unity = y./(repmat(absy,[1,1,size(y,3)])+eps);

res = absy-t;
res = (res + abs(res))/2;
x = unity.*repmat(res,[1,1,size(y,3)]);
