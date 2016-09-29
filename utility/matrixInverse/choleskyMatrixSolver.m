function x = choleskyMatrixSolver(A, b)
% using the cholesky decomposition to solve the linear equation Ax = b
% A: M*N matrix, b : M*K matrix
% first, A is decomposed into L*L' = A, L is a lower triangular matrix
% second, let y = L'*x, solve L*y = b
% finally, solve L'x = y

%% check the matrix consistency
sA = size(A);
sb = size(b);
M = sA(1);
N = sA(2);
K = sb(2);
if ( sA(1) ~= sb(1) )
    error('Matrix size does not match -- A and b');
    x = [];
    return;
end

%% perform the chlesky decomposition
[L, p] = chol(A, 'lower');

if ( p > 0 )
    error('A is not positive definite');
    x = [];
    return;
end

lower.LT = true ;
upper.LT = true ;
upper.TRANSA = true ;
y = linsolve(L, b, lower);
x = linsolve(L, y, upper);
 
% % solve y for L*y = b
% y = zeros(N, K);
% 
% % for every column
% for c=1:K    
%     % starting from the first row, as L is a lower triangler matrix
%     % do the backsubstitution
%     y(1,c) = b(1,c)/L(1,1);    
%     for r=2:N
%         accV = 0;
%         for s=1:r-1
%             accV = accV + y(s,c)*L(r,s);
%         end
%         y(r,c) = (b(r,c)-accV)/L(r,r);
%     end
% end
% 
% % solve x for L'*x = y
% x = zeros(N, K);
% transposeL = L';
% % for every column
% for c=1:K        
%     % starting from the last row
%     % do the backsubstitution
%     x(N,c) = y(N,c)/transposeL(N,N);    
%     for r=N-1:-1:1
%         accV = 0;
%         for s=r+1:N
%             accV = accV + x(s,c)*transposeL(r,s);
%         end
%         x(r,c) = (y(r,c)-accV)/transposeL(r,r);
%     end
% end

