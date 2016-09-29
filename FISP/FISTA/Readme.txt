

fista_slep is for solving the following problem
    min  1/2 \|Ax -y\|^2 + lambda \|Wx\|_1 
It operates as:
         [x, funVal,  ValueL, mse]= fista_slep(A, AT, y, W, WT, lambda, opts);
Here, x is an image (denoted by a vector).





% Function fista_slep_2D
%      FISTA for 2D sum of squares reconstruction
%
%% Problem
%
%  min  1/2 \|Ax -y\|_F^2 + lambda *  \| W x \|_{2,1}
%
%  x=[x(:,1), x(:,2), ..., x(:,c]]
%           x(:,i) is a concatenated vector from the 2D matrix
%           c is the number of coils
%  y=[y_1; y_2; y_c]
%           y_i is a concatenated vector in the undersampled k-space
%  the defulat value for w are all ones
%  x is 2D
%  y is 1D
%  w is 1D
fista_slep_2D is for solving the following problem
    min  1/2 \|Ax -y\|_F^2 + lambda *  \| W x \|_{2,1}
It operates as:
         [x, funVal,  ValueL, mse]=fista_slep_2D(A, AT, y, lambda, opts)
Here, x is a matrix, and corresponds to the image corresponding to different coils.




% Function fista_slep_4D
%      FISTA for 4D sum of squares reconstruction with time
%
%% Problem
%
%  min  1/2 \|Ax -y\|_F^2 + lambda * sum_i \| x(:,:,i).* [w(:,i), w(:,i), ...w(:,i)] \|_{2,1}
%                         + 0.5*rho* sum_i \|x(:,:,i)-x(:,:,i+1)\|_F^2
%
%  x(:,:,i)=[x^i_1, x^i_2, ..., x^i_c]
%           x^i_j is a concatenated vector from the 3D tensor
%           c is the number of coils
%           i goes from 1 to t, with t denoting the number of phases
%  y(index_b:index_e,:)=[y^i_1, y^i_2, ..., y^i_c]
%           y^i_j is a concatenated vector in the undersampled k-space
%           index_b:index_e corresponds to the size of y^i_j
%  the defulat value for w are all ones
%  x is 3D
%  y is 2D
%  w is 2D
It operates as:
    function [x, funVal,  ValueL, mse]=fista_slep_4D(A, AT, y, lambda, opts)
