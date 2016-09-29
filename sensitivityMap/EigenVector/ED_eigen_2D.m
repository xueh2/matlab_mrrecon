function [sensitivityMap, eigD, tau]=ED_eigen_2D(kCalib, kSize, imageSize, opts)
% SVD_eigen_2D:
%    function for computing the Eigen-Vector CSM using the eigen
%    decomposition
%
%% Input parameters: 
%
% kCalib :     calibration data,         [cx, cy, nc)]
% kSize:       calibration kernel size,  [kx, ky]
% imageSize:   the size of the image,    [nx, ny, nc]
% opts:        optional parameters for setting tau
%
%% Output parameters:
%
% sensitivityMap : estimated sensitivity map 
%                 (normalize the output to the phase of the first coil)
% eigD:            the eigenvalues for all the spatial locations
% tau:             the number of kept singular vectors
%
%% Information
%
% Composed by Jun Liu on February 24, 2012
%
% Comments updated on March 7, 2012
%
% Based on Jun Liu's write-up titled
%
% "On the Eigen-Vector Type Approaches 
%   for Coil Sensitivity Maps Estimation"
%
% For any comment/feedback, please contact 
% Jun Liu (jun-liu@siemens.com, junliu.nt@gmail.com)
%
%% Dependency
% This function depends on the function
%    compute_SMH
%


%% extract from opts on how to set tau
% opts.choice =1 ====> 
%              tau is set such that opts.percentage energy is kept
% opts.choice =0 ====>
%              tau is set to the value specified in opts.tau, if provided

if (nargin<4)
    opts=[];
    opts.choice=1;      % set tau to retain opts.percentage energy
    opts.percentage=95; % the percentage of energy to be kept
end

if (~isfield(opts,'choice'))
    opts.choice=1;
end
% if opts.choice is not provided, set opts.choice to 1

if (~opts.choice)
    if (~isfield(opts,'tau'))
        opts.choice=1;
    end
end
% if opts.choice=0, but opts.tau is not provided, set opts.choice=1

if (opts.choice)
    if (~isfield(opts,'percentage'))
        opts.percentage=95;
    end
end
% if opts.choice=1, but opts.percentage is not provided, 
% set opts.percentage=95

%% get the size

nx=imageSize(1);
ny=imageSize(2);
% image size [nx, ny]

nc=size(kCalib,3);
% number of coils

kx=kSize(1);
ky=kSize(2);
% the values should be real and larger than 1

%% generate the matrix A
A=[];
for i=1:nc % 1- nc
    A = [A; im2col(kCalib(:,:,i),[kx ky],'sliding')];
end
% the A matrix
% the elements in one block is stored in the format
% [1  kx+1   ...  ]
% [2  kx+2   ...  ]
% [   ...         ]
% [kx ...     kxky]

%% SVD to obtain v

[u, s, v]=svd(A','econ');
% A'=u * s * v'
s=diag(s);
% get the diagonal entries of s


if (~opts.choice)
    tau=opts.tau;
else
    % set tau such that opts.percentage % energy is retained
    percentage= s./ sum(s) * 100;
    for i=2:length(s)
        percentage(i)=percentage(i-1)+percentage(i);
    end
    tau= min( find(percentage >opts.percentage) );    
end
% set the value for tau


fprintf('\n tau=%d, out of %d \n',tau, min(size(A)));

clear u s A;
v=v(:,1:tau);
% only keep the first tau columns

%% compute SMH 

SMH=compute_SMH(v, kSize, imageSize);
% SMH is a matrix of size [nc, nc, nx, ny]

%% eigen-vector decomposition for each pixel

sensitivityMap=complex(zeros(nc, nx, ny), zeros(nc, nx, ny));
% initialize with zeros
% the reason for using size [nc, nx, ny] is for better memory management in
% Matlab, this matrix is later permuted to [nx, ny, nc]

eigD=zeros(nx,ny);
% initialize eigD with zeros


for j=1:ny
    for i=1:nx
        R = SMH(:, :, i, j);      
        % size: nc x nc
        
        [V, D] = eig( (R + R')/2 );      
        % R is Hermitian theoretically,
        % we apply (R+R')/2 to achieve better numerical stability        
        % eig decomposition
        
        D=diag(D);
        sensitivityMap(:, i, j) = V(:, end);
        eigD(i,j,1)=D(end) / kx/ky;
        eigD(i,j,2)=D(end-1) / kx/ky;
        eigD(i,j,3)=D(end-2) / kx/ky;
        % set the solution as the eigen-vector corresponding to 
        % the largest one 
        %
        % Note that, the eigenvalues are real, which is guaranteed by the
        % Hermitian of R
    end
end

sensitivityMap=permute(sensitivityMap,[2,3,1]);
% permute sensitivityMap to size [nx, ny, nc]

%% normalize the output to the phase of the first coil

aa=conj(sign( sensitivityMap(:,:,1) ));
sensitivityMap=sensitivityMap.* repmat( aa, [1, 1, nc]);
% the csm for the first coil has real and non-negative values

