function [sensitivityMap, eigD, tau]=ED_eigen_3D_parallel(kCalib, kSize, imageSize, opts)
% SVD_eigen_3D:
%    function for computing the Eigen-Vector CSM using the eigen
%    decomposition
%
%% Input parameters:
%
% kCalib :     calibration data,         [cx, cy, cz, nc)]
% kSize:       calibration kernel size,  [kx, ky, kz]
% imageSize:   the size of the image,    [nx, ny, nz, nc]
% opts:        optional parameters for setting tau
%
%% Output parameters:
%
% sensitivityMap : estimated sensitivity map
%                 (normalize the output to the phase of the first coil)
% eigD:            the eigenvalues for all the spatial locations
% tau:             the number of kept singular vectors
%
%% Information:
%
% Composed by Jun Liu on March 2, 2012
%
% Based on Jun Liu's write-up titled
%
% "Revisiting the Eigen-Vector Type Approaches
%   for Coil Sensitivity Maps  Estimation"
%
% For any comment/feedback, please contact
% Jun Liu (jun-liu@siemens.com, junliu.nt@gmail.com)
%
%% Some Notes:
%
% This function yields the same result as ED_eigen_3D (frankly speaking,
% ED_eigen_3D has never been run successfully on my machine, due to
% memory limit; but I have run it conceptually), but it has the
% following key advantages:
%
% 1) The program can be run in parallel (in C++, and maybe Matlab)
% 2) The program requires less memory than ED_eigen_3D
%   (and even less than ED_eigen_2D for the computation memory)
% 3) It is desired that the program can be run on a personal PC
%
%
% The current program assumes that nx >=ny >= nz >=1.
% I believe that, in this setting,
% the program is able to achieve highest computational and storgae
% performance in the C++ enviroment, as the number of fft to be performed
% can be reduced to the minimum.
%
% In other words, if "nx >=ny >= nz" does not satisfy,
% the program can be speeded up by switching the order of nx, ny and nz.
%
%
% This code also include ED_eigen_2D as a special case.
% For this case, please set cz=kz=1 and opts.parallel=1.
%
%
% For computing v from A, SVD is used as a default option. However, when A
% is very large (having both large number of rows and columns), eig on A*A'
% might be a better choice. Of course, from the viewpoint of numerical
% stability, SVD is better. But, I believe that, this is not a big issue,
% especially considering the fact that we only focus on the leading
% eigenvectors/singular vectors.
%
% Other advanced SVD tools can be used for computing v from A, e.g.,
% random PCA.
%
%
% A practical issue with the usage of this approach might be that, how many
% singular vectors to be kept.
%
% Up to now, I use opts.choice=1; and opts.percentage=95;
%    This seems not too back, although I believe that it is not optimal
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
    % 95 might be a too high number
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

%% set the option to get v from A
if (~isfield(opts,'SVD'))
    opts.SVD=1;
end

if (opts.SVD~=0) && (opts.SVD~=1)
    error('\n Currently, only 0 and 1 are supported for opts.SVD');
end
% By default, SVD is used (opts.SVD=1)
%
% However, if A is very large (having both large number of rows and
% columns), eig might be a better choice.
% In this case, please set opts.SVD=0
%
% opts.SVD=2 can be used for random PCA (future work)

%% parallel level
if (~isfield(opts,'parallel'))
    opts.parallel=1;
else
    if (opts.parallel~=1) && (opts.parallel~=2)
        error('\n opts.parallel can only be 1 or 2');
    end
end
% By default, we one level parallel (i.e., in the nz direction)
%
% If memory is a key issue, and two level parallel can be used (in C++)
% In this case, please set opts.parallel=2;


%% get the size

cx=size(kCalib,1);
cy=size(kCalib,2);
cz=size(kCalib,3);
nx=imageSize(1);
ny=imageSize(2);
nz=imageSize(3);
nc=imageSize(4);
kx=kSize(1);
ky=kSize(2);
kz=kSize(3);
% the values should be real and larger than 1
% it is assumed that nx, ny, and nz are all even

if ( (cz==1) && (kz~=1) ) || ( (cz~=1) && (kz==1) )
    error('\n If either cz=1 or kz=1, kz or cz should be 1 too!\n')
end

%% generate the matrix A
%
% fix me, is there any better ways for generating A in the 3D case?
index_i=1:kx;
index_j=1:ky;
index_k=1:kz;
% the index to retrieve the 3D k-space data

A_col=(cx-kx+1)*...
    (cy-ky+1)*...
    (cz-kz+1);
% number of columns in A

A_row=kx*ky*kz* nc;
% number of rows in A

% allocate memory for A
A=single( complex(zeros( A_row, A_col), zeros( A_row, A_col) ) );
%A=complex(zeros( A_row, A_col), zeros( A_row, A_col) );

% generate A
index_put=1;

index_k=1:kz;
for k=1: (cz-kz+1)
    
    index_j=1:ky;
    for j=1:(cy-ky+1)
        
        index_i=1:kx;
        for i=1:(cx-kx+1)
            
            temp=kCalib(index_i,index_j,index_k,:);
            A(:, index_put) = temp(:);
            
            index_put=index_put+1;
            
            index_i=index_i+1;
        end
        index_j=index_j+1;
    end
    index_k=index_k+1;
end
% note that, the order of temp is important!!! in C++ implementation
% is there a better way? (for Matlab implementation)
% for 2D image, im2col can be used
% for C++ implementation, it should be efficient


% the A matrix
% the elements in one block is stored in the format
% [1  kx+1   ...    ]
% [2  kx+2   ...    ]
% [   ...           ]
% [kx ...     kxkykz]

%% SVD to obtain v

fprintf('\n SVD decomposition for v... ');

tic;

if opts.SVD==1
    [u, s, v]=svd(A','econ');
    % A'=u * s * v'
    s=diag(s);
    % get the diagonal entries of s
end

if opts.SVD==0
    [v, s]=eig(A*A');
    % A*A'= v* s.^2 * v'
    
    v=v(:,end:-1:(end-min(size(A))+1));
    % only take the leading eigenvectors
    
    s=sqrt( abs(diag(s) ));
    % it is real theoretically, but the solver may return complex value
    
    s=s(end:-1:(end-min(size(A))+1));
    % only take the eigenvalues
end
% we provide two choices:
%  1) SVD
%  2) eig decomposition
%
% In the future, we would like to make use of random PCA, which seems to be
% faster. Of courese, a question is how to set tau?
%
% Some test need to be done in this direction.
% For 2D, 150 is enough?
% For 3D, 500 is enough?
%
% The percentage of 95% seems to lead to ~200 for the 2D case
%                                        ~1500 for the 3D case
%
% Note that, in theory, if tau is too small, 
%            the performance cannot be good (too smooth)
%

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

toc;

%% compute SMH in a parallel way
% conceptually, SMH is a matrix of size [nc, nc, nx, ny, nz]
% we did not store it explicitly, as it is memory consuming

P=v*v';
% the P matrix
clear v;
% save memory

kkx=kx*2-1;
kky=ky*2-1;
kkz=kz*2-1;
% this is the size of kernel in computing SMH

%% creat the index for generating kernel_nz
B=zeros(kkx,kky,nz);
B(1:kx,1:ky,1:kz)=1;
index_base=find(B);
% a kx* ky* kz vector
% this is used to construct kernel_nz

index_all=zeros(kx*ky*kz, kx*ky*kz);
for i=1:kx
    for j=1:ky
        for k=1:kz
            
            index_all(:,i + (j-1) * kx + (k-1) * kx * ky)=...
                index_base+ ...
                kx -i + ...
                (ky -j)* kkx + ...
                (ceil((nz+1)/2) -k) * kkx * kky;
            
        end
    end
end
% index_all contains the indices to
% for putting the columns of P in appropriate positions
%
% For the 3D image of size [nx,ny,nz],
% the k-space center is [ceil((nx+1)/2), ceil((ny+1)/2),  ceil((nz+1)/2)].
%
% For the 3D image of size [kkx,kky,nz],
% the k-space center is [kx, ky,  ceil((nz+1)/2)].
%
% It is assumed that nx and ny are both even.

% % nx=256; ny=nx; kx=5; ky=5; kz=2;nz=6;
% % % this is used for test the correctness of the code
% % % put the first line before the initialization of A
% %
% % B=zeros(kkx,kky,nz); B(index_all(:,1))=1;
% % for i=1:(kx*ky*kz) B(index_all(:,i))=1; end
% % % test the correctness of index_all

%% generate kernel_nz

% as MSH (see my write-up) is Hermitian, we only need to store the
% upper-triangular part of kernel_nz

fprintf('\n Generating the decoupled kernel along nz direction... ');

tic;

ncnc=nc*(nc+1)/2;
% ncnc is about half of nc*nc
% this shall save the memory cost and computation cost by roughly half

kernel_nz=complex(zeros(kkx,kky,nz, ncnc), zeros(kkx,kky,nz, ncnc));
% kernel_nz decouples nz direction
% of size [kkx,kky,nz, ncnc]

index_put=1;
% the index to put entries in the ncnc direction of kernel_nz(:,:,:,ncnc)

for i=1:nc % i-th column, see Section 4.1
    for j=1:i %j-th row, see Section 4.1
        
        kernel=complex(zeros(kkx, kky, nz), zeros(kkx, kky, nz));
        % this is kernel to filled in kernel_nz
        
        index_j=(1:(kx*ky*kz))+ (j-1)*kx*ky*kz;
        % the index of P which contributes to the j-th row of SMH
        
        for t=1:(kx*ky*kz)
            kernel(index_all(:,t))=kernel(index_all(:,t)) +...
                P(index_j, t+ (i-1)*kx*ky*kz );
        end
        % fill in the summation of the kernel
        
        kernel_nz(:,:,:,index_put)=kernel;
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
        % conceptually, kernel_nz is of size [kkx,kky,nz,nc,nc]
        % we only fill in the upper-triangular part.
        % we leave lower-triangular of kernel_nz unfilled, as this part can
        % be obtained from its upper-triangular part
    end
end

kernel_nz=fftshift(...
    ifft(...
    fftshift(kernel_nz,3),...
    [],3), 3) * sqrt(nz);
% kernel_nz is of size [kkx,kky,nz, ncnc]
% apply ifft along nz direction
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

toc;

%% initialize memory for our results in sensitivityMap and eigD
sensitivityMap=complex(zeros(nc, nx, ny, nz), zeros(nc, nx, ny, nz));
% initialize with zeros
% the reason for using size [nc, nx, ny, nz] is for better memory management in
% Matlab, this matrix is later permuted to [nx, ny, nc]

eigD=zeros(nx,ny,nz);
% initialize eigD with zeros


%% Apply a parallel implementation along the nz direction
if opts.parallel==1
    for i_nz= 1:nz
        
        fprintf('\n parallel in nz direction: %d out of %d... ',i_nz,nz);
        tic;
        
        kernel=complex(zeros(nx, ny, ncnc), zeros(nx, ny, ncnc));
        % of size [nx, ny, ncnc]
        
        kernel( (ceil((nx+1)/2) - (kx-1)): (ceil((nx+1)/2) + (kx-1)),...
            (ceil((ny+1)/2) - (ky-1)): (ceil((ny+1)/2) + (ky-1)),:)=...
            squeeze(kernel_nz(:,:,i_nz,:));
        % for nx and ny direction, put kernel_nz(:,:,i_nz,:,:) at the center
        %
        % as can seen above, we do not decouple nx or ny direction,
        % to further save memory and/or computation cost (in C++), please
        % go the setting opts.parallel=2
        
        
        if (mod(ncnc,2)==0)
            SMH=fftshift(...
                ifft2(fftshift(kernel) ) )*...
                sqrt(nx*ny);
        else
            SMH=fftshift( fftshift(...
                ifft2(fftshift(fftshift(kernel,1),2)),...
                1),2 )*...
                sqrt(nx*ny);
        end
        % SMH is of size [nx, ny, ncnc]
        % ifft2 to the image space along the third direction (ncnc)
        %
        % note that, the "else" case might be the one used in C++
        % implementation, the "if" case takes advantage of the fast
        % computation of ifft2 along the "ncnc" direction, and that fftshift
        % and fftshift along the ncnc direction can be canceled due to that
        % ncnc is even
        
        SMH=permute(SMH,[3,1,2]);
        % permute to [ncnc, nx, ny]
        
        sensitivityMap_inz= complex(zeros(nc, nx, ny), zeros(nc, nx, ny));
        % the csm obtained along the nz direction
        
        for j=1:ny
            for i=1:nx
                
                R=complex( zeros(nc,nc), zeros(nc,nc) );
                % set R to be a complex zero matrix of size [nc,nc]
                
                R(index_R) = SMH(:, i, j);
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
                % note the following structure of SMH, kernel and kernel_nz
                % [1 2 4  ...     ]
                % [  3 5  ...     ]
                % [    6  ...     ]
                % [          ncnc ]
                
                temp= R + triu(R,1)';
                % fill in the low-triangular of R
                
                [V, D] = eig( (temp + temp')/2 );
                % to make sure that the input to eig is Hermitian
                % In C++, (temp + temp')/2 can be combined with filling R
                % with SMH
                %
                % theoretically, temp should be Hermitian (the diagonal
                % entries of R and temp should be real
                %
                % but minor numerical error might lead to significantly
                % different solution
                
                sensitivityMap_inz(:, i, j) = V(:, end);
                eigD(i,j,i_nz)=D(end)/ kx/ky/kz;
                % scale by kx*ky*kz, See Section
            end
            
            
            sensitivityMap(:,:,:,i_nz)=sensitivityMap_inz;
            % copy sensitivityMap_inz to the matrix sensitivityMap
        end
        toc;
    end
    % the nz direction is decoupled
    % parallel implementation for the loop i_nz
end

%% Apply a parallel implementation along the nz and ny direction
if opts.parallel==2
    
    for i_nz= 1:nz
        
        fprintf('\n parallel in nz direction: %d out of %d... ',i_nz,nz);
        tic;
        
        %% generate the kernel_ny from kernel_nz (size [kkx,kky,i_nz, ncnc])
        
        kernel_ny=complex(zeros(kkx, ny, ncnc), zeros(kkx, ny, ncnc));
        % of size [kkx, ny, ncnc]
        
        kernel_ny(:,...
            (ceil((ny+1)/2) - (ky-1)): (ceil((ny+1)/2) + (ky-1)),...
            :)=...
            squeeze(kernel_nz(:,:,i_nz,:));
        % for ny direction, put kernel_nz(:,:,i_nz,:) at the center
        %
        % as can seen above, we decouple ny direction,
        
        kernel_ny=fftshift(...
            ifft(...
            fftshift(kernel_ny,2),...
            [],2), 2) * sqrt(ny);
        % kernel_ny is of size [kkx,ny, ncnc]
        % apply ifft along ny direction
        % this direction now has been decoupled
        
        
        sensitivityMap_inz= complex(zeros(nc, nx, ny), zeros(nc, nx, ny));
        % the csm obtained along the nz direction
        
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
            
            SMH_ny=permute(SMH_ny,[2,1]);
            % permute to [ncnc, nx]
            
            sensitivityMap_iny=complex(zeros(nc, nx), zeros(nc, nx));
            % the csm obtained along the ny direction (for a given inz)
            
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
                % note the following structure of SMH, kernel and kernel_nz
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
                eigD(i,i_ny,i_nz)=D(end)/ kx/ky/kz;
                % scale by kx*ky*kz, See Section
            end
            
            sensitivityMap_inz(:, :, i_ny)=sensitivityMap_iny;
            % sensitivityMap_iny of size [nc, nx]
            % sensitivityMap_inz of size [nc, nx, ny]
        end
        % parallel along ny ends here
                
        sensitivityMap(:,:,:,i_nz)=sensitivityMap_inz;
        % copy sensitivityMap_inz to the matrix sensitivityMap
        % sensitivityMap of size [nc, nx, ny, nz]
        toc;
    end
    % the nz direction is decoupled
    % parallel implementation for the loop i_nz ends here
end



%% normalize the output to the phase of the first coil

sensitivityMap=permute(sensitivityMap,[2, 3, 4, 1]);
% permute sensitivityMap to size [nx, ny, nz, nc]

aa=conj(sign( sensitivityMap(:,:,:,1) ));
sensitivityMap=sensitivityMap.* repmat( aa, [1, 1, 1, nc]);
% the csm for the first coil has real and non-negative values

