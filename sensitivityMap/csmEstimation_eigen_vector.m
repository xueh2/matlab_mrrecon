function [sensitivityMap, eigValue, tau]=csmEstimation_eigen_vector(k_center, imageSize, kernelSize, tau)
%% Estimate csm using the eigen vector approach
%
% Composed by Jun Liu on January 20, 2012
% For any problem/question, please contact Jun Liu (jun-liu@siemens.com)
%
% Input: 
%   k_center  -- the k-space data in the center
%   imageSize  --the size of image [fe, pe, c (channel)]
%   kernelSize --the size of the GRAPPA kernel (block size)
%               - we assume that kernelSize is even 
%               - fix me, for odd number
%

%% generate the matrix A
B=[];
for i=1:imageSize(3)
    B = [B; im2col(k_center(:,:,i),[kernelSize kernelSize],'sliding')];
end

%% SVD to obtain v
A=B';
[u, s, v]=svd(A,'econ');
s=diag(s);

% set tau such that 90% energy is retained
percentage= s./ sum(s) * 100;
for i=2:length(s)
    percentage(i)=percentage(i-1)+percentage(i);
end
tau= min( find(percentage >90) );
% either tau or 90% may be used as a parameter later
%
tau
%sum(s(1:tau))/sum(s)


%% convert to image space
x_start=imageSize(1)/2 -kernelSize/2+1;
x_end=imageSize(1)/2 + kernelSize/2;
y_start=imageSize(2)/2 -kernelSize/2+1;
y_end=imageSize(2)/2 + kernelSize/2;

index=1:kernelSize^2;
for i=1:imageSize(3)
    for j=1:tau
        aa=reshape( v(index,j), [kernelSize,kernelSize]);
        
        AA=zeros([imageSize(1), imageSize(2)]);
        AA(x_start: x_end, y_start: y_end)=aa;
        
        imageAll(:,:, i, j)=fftshift( ifft2(fftshift(AA)) )* sqrt( prod(imageSize(1:2)) );
    end
    
    index=index+kernelSize^2;    
end


%  figure; imshow( abs(imageAll(:,:,1,1)),[])
%  figure; imshow( abs(imageAll(:,:,1,2)),[])
%  figure; imshow( abs(imageAll(:,:,1,3)),[])
%  
%  figure; imshow( abs(imageAll(:,:,2,1)),[])
%  figure; imshow( abs(imageAll(:,:,2,2)),[])
%  figure; imshow( abs(imageAll(:,:,2,3)),[])

% imageAll is of size [fe,pe, tau, channel]

%% eigen-vector decomposition for each pixel
for i=1:imageSize(1)
    for j=1:imageSize(2)
        R = squeeze(imageAll(i, j, :, :)); % size: c x tau
        [V, D] = eig( R*R' );              % size: c x c
        sensitivityMap(i, j, :) = V(:, end);
        eigValue(i,j)=D(end);
    end
end

%% rescale by a scalar to make sure that the coil sensitivity map for the last coil has positive (real) values
aa=sign( sensitivityMap(:,:,end) );

sensitivityMap=sensitivityMap.* repmat( aa, [1, 1, imageSize(3)]);
