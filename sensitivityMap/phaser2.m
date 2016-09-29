clf

% Load the image domain
%load ../data/human.mat
%load ../data/phantom.mat
% Get the sizes
[Nx,Ny,Nc] = size(data);

% % Whiten the data
% N = reshape(data(1:50,1:50,:),[50*50 Nc]);
% Cni = (1/(size(N,1)-1)) * (N'*N);
% L = chol(Cni,'lower');
% Linv = inv(L);
% data = reshape((reshape(data,[Nx*Ny Nc]) * Linv),[Nx Ny Nc]);
% % At this point, the noise variance in each channel is 1
% % and we can use this for regularization.
% N = reshape(data(1:50,1:50,:),[50*50 Nc]);
% Cnw = (1/(size(N,1)-1)) * (N'*N);

if (0)
    % create a random unitary matrix
    H = randn(Nc)+1j*randn(Nc);
    [U,S,V] = svd(H);
    data = reshape((reshape(data,[Nx*Ny Nc]) * V),[Nx Ny Nc]);    
    N = reshape(data(1:50,1:50,:),[50*50 Nc]);
    Cnv = (1/(size(N,1)-1)) * (N'*N);
end

if (0)
    D = reshape(data,[Nx*Ny Nc]);
    [U,S,V] = svd(D,'econ');
    datav = reshape(U*S,[Nx Ny Nc]);
    phiref = angle(squeeze(sum(sum(datav,1),2)));
    % for p = 1:32
    %     datac(:,:,p) = exp(-1j*phiref(p))*datav(:,:,p);
    % end
    N = reshape(datav(1:50,1:50,:),[50*50 Nc]);
    Cns = (1/(size(N,1)-1)) * (N'*N);
    data = datav;
end

if (0)
    % Phase each channel wrt its mean
    phiref = angle(squeeze(sum(sum(data,1),2)));
    for p = 1:32
        data(:,:,p) = exp(-1j*phiref(p))*data(:,:,p);
    end
    N = reshape(data(1:50,1:50,:),[50*50 Nc]);
    Cnm = (1/(size(N,1)-1)) * (N'*N);
end
    
% Smoothness scale
scale = 7;

% The number of iterations
Niter = 10;

% Type of phase constraint
% 0 global phase
% 1 high pass
% 2 sum(B) is real
% otherwisef none
phaseConstraint = 0;

% Initialization
starttype = 3;
switch starttype
    case 0
        % Initial guess for R is SOS
        R = sqrt(sum(abs(data).^2,3));

    case 1
        % Initial guess is sum of the data
        % this is generally a bad idea
        R = sum(data,3);

    case 2
        % Initial guess is random numbers
        R = randn(Nx,Ny) + 1j*randn(Nx,Ny);

    otherwise
        % Assume constant coil sensitivities
        b = squeeze(sum(sum(data,1),2));
        b = b/norm(b);
        R = zeros(Nx,Ny);
        for c = 1:Nc
            R = R + conj(b(c))*data(:,:,c);
        end
end

% Store the initial guess
Ri = R;

% Display some results
figure(1)
subplot(1,4,1); imagesc(abs(R)); axis image; colorbar
title(sprintf('Abs(R) scale: %d, iter: %d',scale,0))
subplot(1,4,2); imagesc(angle(R)); axis image; colorbar
title(sprintf('Angle(R) scale: %d, iter: %d',scale,0))
subplot(1,4,3); imagesc(zeros(Nx,Ny),[-pi pi]); axis image; colorbar
title(sprintf('Phi scale: %d, iter: %d',scale,0))
subplot(1,4,4); imagesc(zeros(Nx,Ny)); axis image; colorbar
title(sprintf('Abs(B) scale: %d, iter: %d',scale,0))
drawnow;

%pause(1.0)
pause

% The convolution kernels
G = ones(scale); G = G/sum(G(:));

for iter = 1:Niter

    % v1*s1^2 i.e. (u1*s1)'*D
    B = repmat(conj(R),[1 1 Nc]).*data;
    for c = 1:Nc
        B(:,:,c) = conv2(B(:,:,c),G,'same');
    end
    
    % v1 i.e.
    N = sqrt(sum(abs(B).^2,3));
    B = B ./ repmat(N,[1 1 Nc]);

    % u1*s1, i.e. D*v1'
    R = sum(conj(B).*data,3);
    
    switch phaseConstraint
        case 0
            % Weighted sum of all the coil sensitivities
            dhat = repmat(R,[1 1 Nc]).*B;
            b = squeeze(sum(sum(dhat,1),2));
            b = b/norm(b);
            T = zeros(Nx,Ny);
            for c = 1:Nc
                T = T + conj(b(c))*B(:,:,c);
            end
            T = angle(T); 
            % this is all the phase that is not already in R
            
            % Put all of the phase variation into the image
            R = exp(1j*T).*R;
            for c = 1:Nc
                B(:,:,c) = exp(-1j*T).*B(:,:,c);
            end

        case 1
            Rc = conv2(R,G,'same');
            T = -angle(Rc);
            R = exp(1j*T).*R;
            for c = 1:Nc
                B(:,:,c) = exp(-1j*T).*B(:,:,c);
            end
            
    end

    % This is NOT the same as
    % Rc = sum(conj(Bc).*data,3);    
    % But it is if you  phase the data as well
    % this is phase each channel wrt to it's mean
    % value and say that the phase of the sum
    % of the coil sensitivities is zero.
    % We have to include R because it also has some
    % smooth phase that is not in B and we have to take
    % that into account.
    
    % Display some results
    figure(1)
    subplot(1,4,1); imagesc(abs(R)); axis image; colorbar
    title(sprintf('Abs(R) scale: %d, iter: %d',scale,iter))
    subplot(1,4,2); imagesc(angle(R)); axis image; colorbar
    title(sprintf('Angle(R) scale: %d, iter: %d',scale,iter))
    subplot(1,4,3); imagesc(T,[-pi pi]); axis image; colorbar
    title(sprintf('Phi scale: %d, iter: %d',scale,iter))
    subplot(1,4,4); imagesc(sum(abs(B),3)); axis image; colorbar
    title(sprintf('Abs(B) scale: %d, iter: %d',scale,iter))
    drawnow;
    
    %pause(1.0);

end

% Recenter the coil sensitivities
% the object indepent sum
h = squeeze(sum(sum(B,1),2));
h = h/norm(h);
T = zeros(Nx,Ny);
for c = 1:Nc
    T = T + conj(h(c))*B(:,:,c);
end
T = angle(T);
Bf = B;
for c = 1:Nc
    Bf(:,:,c) = exp(-1j*T).*B(:,:,c);
end
Rf = sum(conj(Bf).*data,3);
