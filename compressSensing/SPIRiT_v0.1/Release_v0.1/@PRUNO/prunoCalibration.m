function kernel = prunoCalibration(kCalib,kSize,thresEigenValue, r, show)
% compute the pruno kernel
% kspace: the fullkspace (Nfe Npe Coil), calibration area
% kSize : kernel size, kSize(1) is along row direction, kSize(2) along column
% thresEigenValue : threshold to get the nullspace, relative to the maximal eigenvalue; if <0, r null space kernels are selected
% show : if 1, plot a figure to show the eigenvalue
% kernel : a kSize(1)*kSize(2)*nCoil*r matrix, corresponding to r pruno kernels

disp(['calibration region size : ' num2str(size(kCalib))]);
disp(['kernel size : ' num2str(kSize)]);

[sx,sy,nCoil] = size(kCalib);

A = [];

dummyK = zeros(kSize(1),kSize(2)); 
ind = find(dummyK==0);

% composite the data matrix
for n=1:nCoil
	tmp  = im2col(kCalib(:,:,n),kSize,'sliding').';
	A = [A; tmp'];
end
AAt = A*A';

% perform SVD
[U,S,V] = svd(AAt);

N = kSize(1)*kSize(2)*nCoil;
nullSpaceCutOff = r; % the number of null space vectors
E = sqrt(abs(diag(S))); % the eigenvalues are squared after A*A'

if ( thresEigenValue > 0 )   
    v = thresEigenValue*E(1);
    ind = find(E<v);
    if ( ~isempty(ind) ) 
        nullSpaceCutOff = N - ind(1);
    else
        nullSpaceCutOff = floor(N/3);
    end
end

% get the null space kernel
kernel = U(:, end-nullSpaceCutOff+1:end);
kernel = reshape(kernel, [kSize(1) kSize(2) nCoil nullSpaceCutOff]);

% show results
if ( show )    
    figure 
    hold on
    plot(E);
    plot(E, 'rx');
    plot(N-nullSpaceCutOff+1:N, E(end-nullSpaceCutOff+1:end), 'ko');
    hold off
       
    figure; semilogy(E);
    hold on
    semilogy(E, 'rx');
    semilogy(N-nullSpaceCutOff+1:N, E(end-nullSpaceCutOff+1:end), 'ko');
    hold off
end


