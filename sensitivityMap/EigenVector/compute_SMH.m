function  SMH=compute_SMH(v, kSize, imageSize)
% compute_SMH:
%    function for computing S * M'
%
% For each spatial location, M is of size [nc, nc*kx*ky]
%   [ 1 1 ... 1                               ]
%   [            1 1 ..., 1                   ]
%   [ ...                                     ]
%   [                        ...   1 1 ..., 1 ]
%
% For each spatial location, S is also of size [nc, nc*kx*ky]
%
%% Input parameters: 
%
% v :          singular vectors,         [nc*kx*ky, tau]
% kSize:       calibration kernel size,  [kx, ky]
% imageSize:   the size of the image
%
%% Output parameters:
%
% SMH:         S * M' (size [nc, nc, nx, ny])
%
%% Information
%
% Composed by Jun Liu on February 24, 2012
%
% Based on Jun Liu's write-up titled
%
% "Revisiting the Eigen-Vector Type Approaches 
%   for Coil Sensitivity Maps  Estimation"
%
% For any comment/feedback, please contact Jun Liu 
% jun-liu@siemens.com, junliu.nt@gmail.com
%
% As SMH is Hermitian (for each pixel loation),
% the computation cost can be saved by half.
% See Theorem 1 of the write-up
% Updated on Feburary 27, 2012


%% get the size

nx=imageSize(1);
ny=imageSize(2);
nc=imageSize(3);
kx=kSize(1);
ky=kSize(2);
% the values should be real and larger than 1
% checking is not conducted here

P=v*v';
% the P matrix

%% creat an index array for placing the values of in proper positions

A=zeros(nx,ny);
A(1:kx,1:ky)=1;
index_base=find(A);
% a kx* ky vector

index_all=zeros(kx*ky, kx*ky);
for i=1:kx
    for j=1:ky
        index_all(:,i + (j-1) * kx)=index_base+ ...
            ceil((nx+1)/2) -i + ...
            (ceil((ny+1)/2) -j) * nx;
    end
end
% index_all contains the indices to
% for putting the columns of P in appropriate positions
%
% For the 2D image of size [nx,ny],
% the k-space center is [ceil((nx+1)/2), ceil((ny+1)/2)].
%
% It is assumed that nx and ny are both even.

% % nx=256; ny=nx; kx=5; ky=5;
% % % this is used for test the correctness of the code
% % % put the first line before the initialization of A
% % 
% % B=zeros(nx,ny); B(index_all(:,25))=1;
% % % test the correctness of index_all


%% fill in SMH

SMH=complex(zeros(nx, ny, nc, nc), zeros(nx, ny, nc, nc));
% initialize SMH with zeros
% for better memory management in Matlab, 
% [nx,ny,nc,nc] is used in computation

for i=1:nc % i-th column, see Section 4.1
    for j=1:i %j-th row, see Section 4.1
        
        kernel=complex(zeros(nx, ny), zeros(nx, ny));
        % this is kernel to be applied ifft2 later
        
        index_j=(1:(kx*ky))+ (j-1)*kx*ky;
        % the index of P which contributes to the j-th row of SMH
        
        for t=1:(kx*ky)
            kernel(index_all(:,t))=kernel(index_all(:,t)) +...
                P(index_j, t+ (i-1)*kx*ky );
        end
        % fill in the summation of the kernel
        
        SMH(:,:,j,i)=fftshift(ifft2(fftshift(kernel)))* sqrt(nx*ny);
        % fill the j-th row and i-th column of SMH
        % see write-up (Section 4.1)
        
        if (j~=i)
            SMH(:,:,i,j)=conj(SMH(:,:,j,i));
            % make use of the Hermitan property of SMH
        end
    end
end

SMH=permute(SMH,[3,4,1,2]) * sqrt(nx*ny);
% permute to [nc, nc, nx, ny]
% "*sqrt(nx*ny)" is added
% the reason to multiply by "sqrt(nx*ny)" is due to the covolution
% operator in Section 2.2

