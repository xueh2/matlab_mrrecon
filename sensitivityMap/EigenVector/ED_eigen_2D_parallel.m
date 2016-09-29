function [sensitivityMap, eigD, tau]=ED_eigen_2D_parallel(kCalib, kSize, imageSize, opts)
% ED_eigen_2D_parallel:
%    function for computing the Eigen-Vector CSM using the eigen
%    decomposition (in a parallel mode along ny direction)
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
% Composed by Jun Liu on March 7, 2012
%
% Based on Jun Liu's write-up titled
%
% "Revisiting the Eigen-Vector Type Approaches 
%   for Coil Sensitivity Maps  Estimation"
%
% For any comment/feedback, please contact 
% Jun Liu (jun-liu@siemens.com, junliu.nt@gmail.com)
%
%% Notes
% This function should yield the same solution as ED_eigen_2D


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

nc=size(kCalib,3);

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

%% compute SMH in a parallel way (preparation stage)
% conceptually, SMH is a matrix of size [nc, nc, nx, ny]
% we did not store it explicitly, as it is memory consuming

P=v*v';
% the P matrix
clear v;
% save memory

kkx=kx*2-1;
% kkx is the size of kernel in the nx direction

%% creat the index for generating kernel_ny
B=zeros(kkx,ny);
B(1:kx,1:ky)=1;
index_base=find(B);
% a kx* ky vector

index_all=zeros(kx*ky, kx*ky);
for i=1:kx
    for j=1:ky
        index_all(:,i + (j-1) * kx)=index_base+ ...
            kx -i + ...
            (ceil((ny+1)/2) -j) * kkx;
    end
end
% index_all contains the indices to
% for putting the columns of P in appropriate positions


%% generate kernel_ny

% as MSH (see my write-up) is Hermitian, we only need to store the
% upper-triangular part of kernel_ny

%fprintf('\n Generating the decoupled kernel along ny direction... ');

ncnc=nc*(nc+1)/2;
% ncnc is about half of nc*nc
% this shall save the memory cost and computation cost by roughly half

kernel_ny=complex(zeros(kkx,ny, ncnc), zeros(kkx,ny, ncnc));
% kernel_ny decouples ny direction
% of size [kkx,ny, ncnc]

index_put=1;
% the index to put entries in the ncnc direction of kernel_ny(:,:,:,ncnc)

for i=1:nc % i-th column, see Section 4.1
    for j=1:i %j-th row, see Section 4.1
        
        kernel=complex(zeros(kkx, ny), zeros(kkx, ny));
        % this is kernel to filled in kernel_ny
        
        index_j=(1:(kx*ky))+ (j-1)*kx*ky;
        % the index of P which contributes to the j-th row of SMH
        
        for t=1:(kx*ky)
            kernel(index_all(:,t))=kernel(index_all(:,t)) +...
                P(index_j, t+ (i-1)*kx*ky );
        end
        % fill in the summation of the kernel
        
        kernel_ny(:,:,index_put)=kernel;
        % kernel_ny is of size [kkx,ny, ncnc]
        % kernel    is of size [kkx, ny]
        % i-th column, j-th row, 
        % stored in the following way
        % [1 2 4  ...     ]
        % [  3 5  ...     ]
        % [    6  ...     ]
        % [          ncnc ]
        %
        % this order is important for the eigen-decomposition later
        
        index_put=index_put+1;
        % increment the index by one
        
        % note:
        % conceptually, kernel_ny is of size [kkx,ny,nc,nc]
        % we only fill in the upper-triangular part.
        % we leave lower-triangular of kernel_ny unfilled, as this part can
        % be obtained from its upper-triangular part
    end
end

kernel_ny=fftshift(...
    ifft(...
    fftshift(kernel_ny,2),...
    [],2), 2) * sqrt(ny);
% kernel_ny is of size [kkx,ny, ncnc]
% apply ifft along ny direction
% this direction now has been decoupled
%

C=triu(ones(nc,nc));
% C is
% [1 1 1  ...     ]
% [  1 1  ...     ]
% [    1  ...     ]
% [             1 ]

index_R=find(C);
% index_R contains the index to fill in the matix (see later)
%
% [1 2 4  ...     ]
% [  3 5  ...     ]
% [    6  ...     ]
% [          ncnc ]


%% Apply a parallel implementation along the ny direction

for i_ny=1:ny
    %fprintf('\n parallel in ny direction: %d out of %d... ',i_ny,ny);
    
    kernel=complex(zeros(nx, ncnc), zeros(nx, ncnc));
    % of size [nx, ncnc]
    
    kernel( (ceil((nx+1)/2) - (kx-1)): (ceil((nx+1)/2) + (kx-1)),...
        :)=...
        squeeze(kernel_ny(:,i_ny,:));
    % for nx direction, put kernel_ny(:,i_ny,:) at the center
    
    SMH_ny=fftshift(...
        ifft(fftshift(kernel,1),[],1),1 )*...
        sqrt(nx);
    % SMH is of size [nx, ncnc]
    % ifft to the image space along the second direction (ncnc)
    
    SMH_ny=permute(SMH_ny,[2,1])* sqrt(nx*ny);
    % permute to [ncnc, nx]
    % the reason to multiply by "sqrt(nx*ny)" is due to the covolution
    % operator in Section 2.2
    
    sensitivityMap_iny=complex(zeros(nc, nx), zeros(nc, nx));
    % the csm obtained along the ny direction
    
    for i=1:nx
        R=complex( zeros(nc,nc), zeros(nc,nc) );
        % set R to be a complex zero matrix of size [nc,nc]
        
        R(index_R) = SMH_ny(:, i);
        % size: ncnc
        % put the entries of SMH to the upper-triangular of R
        %
        % index_R is the indices for
        %
        % [1 1 1  ...     ]
        % [  1 1  ...     ]
        % [    1  ...     ]
        % [             1 ]
        %
        % note the following structure of SMH, kernel and kernel_ny
        % [1 2 4  ...     ]
        % [  3 5  ...     ]
        % [    6  ...     ]
        % [          ncnc ]
        
        temp= R + triu(R,1)';
        % fill in the low-triangular of R
        
        [V, D] = eig( (temp + temp')/2 );
        % to make sure that the input to eig is Hermitian
        %
        % theoretically, temp should be Hermitian (the diagonal
        % entries of R and temp should be real
        %
        % but minor numerical error might lead to significantly
        % different solution
        
        sensitivityMap_iny(:, i) = V(:, end);
        eigD(i,i_ny)=D(end) / kx/ky;
        % scale by kx*ky, See Section
    end
    
    sensitivityMap(:, :, i_ny)=sensitivityMap_iny;
    % sensitivityMap_iny of size [nc, nx]
    % sensitivityMap of size [nc, nx, ny]
end
% parallel along ny ends here

sensitivityMap=permute(sensitivityMap,[2,3,1]);
% permute sensitivityMap to size [nx, ny, nc]

%% normalize the output to the phase of the first coil

aa=conj(sign( sensitivityMap(:,:,1) ));
sensitivityMap=sensitivityMap.* repmat( aa, [1, 1, nc]);
% the csm for the first coil has real and non-negative values

