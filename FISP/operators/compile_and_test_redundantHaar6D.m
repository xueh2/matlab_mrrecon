mex -setup
%choose your best compiler in the list
%if you see that the openmp optimization is not working (i.e. not 100% of the
%processors are used), change the options file
%(...\Application Data\MathWorks\MATLAB\R2008a\mexopts.bat) to add the
%COMPFLAG /openmp

mex redundantHaar6D.c COMPFLAGS="$COMPFLAGS /openmp"
mex redundantHaarPseudoProx.c COMPFLAGS="$COMPFLAGS /openmp"
mex wavL1norm.c COMPFLAGS="$COMPFLAGS /openmp"
mex wavSoftTreshold.c  COMPFLAGS="$COMPFLAGS /openmp"


%% Data
cd ../core
load_libraries

x=(rand(128,128,1,1,24)+1i*rand(128,128,1,1,24));
%x=(rand(128,128,1,1,24); %to test is all files are correct, you may need
%to perform the tests for complex matrices and for real matrices.
W=redundantHaarOperator(1)*redundantHaarOperator(2)*redundantHaarOperator(5);
Wx=W*x;
al=@(x) x(:);

%% first test, C++ Vs Matalab for wavelets
%%
%create data and appropriate matlab haar operator




% test forward
fprintf('Norm of the difference between the 2 forward implementations: %f.\n',norm(al(W*x-redundantHaar6D(x))))
tic
for k=1:10
    W*x;
end
dt=toc();
fprintf('Matlab implementation performs forward wavelets in %fs.\n',dt)
tic
for k=1:10
    redundantHaar6D(x);
end
dt=toc();
fprintf('C++ implementation performs forward wavelets in %fs.\n',dt)

% test backward
fprintf('\nNorm of the difference between the 2 backward implementations: %f.\n',norm(al(W'*Wx-redundantHaar6D(Wx,-1))))
tic
for k=1:10
    W'*Wx;
end
dt=toc();
fprintf('Matlab implementation performs backward wavelets in %fs.\n',dt)
tic
for k=1:10
    redundantHaar6D(Wx,-1);
end
dt=toc();
fprintf('C++ implementation performs backward wavelets in %fs.\n',dt)


%% Proximal operator and L1 norm implementation tests
%% 
fprintf('\n')

%use C++ implementation to go faster
W=SimpleOperator(@(x) redundantHaar6D(x), @(x) redundantHaar6D(x,-1));

%weights
time_factor=2;
lambda=0.1;

%create the weights for C++ implementation
lambdas=lambda*ones(2,2,2,2,2,2);
lambdas(:,:,:,:,2,:)=time_factor*lambdas(:,:,:,:,2,:);

%create weighted (pseudo) proximal operator for matlab implementation
%(be careful, in pure theory, as the wavelets are redundant, this is NOT
%the exact proximal operator)
S=scaleHaarWav(size(W*zeros(size(x))),5,time_factor);
prox=@(x,mu) W'*S*(softThresholding(S'*W*x,mu));

% BEWARE:
% Matlab and C++ implementation cannot be compared anymore,
%indeed, the last element in high frequency is no more tresholded in the
%C++ implementation as it is the difference of 1st and last element and
%then does not have any interest. Its cost is left to 0.



%% L1 norm, proximal and tresholding tests

fprintf(1,'Difference between using wavelets tresholding and reprojecting,\n and using proximal implementation: %e.\n',norm(al(W'*wavSoftTreshold(W*x,lambdas)-redundantHaarPseudoProx(x,lambdas))))

x=(rand(3,4,3)+1i*rand(3,4,3));
tic
for k=1:5
    redundantHaarPseudoProx(x,lambdas);
end
dt=toc();
fprintf('Proximal implementation needs %fs.\n',dt)
tic
for k=1:5
    W'*wavSoftTreshold(W*x,lambdas);
end
dt=toc();
fprintf('Wavelet tresholding and projecting implementation needs %fs.\n',dt)


[prox_x,L1norm]=redundantHaarPseudoProx(x,lambdas);

fprintf('Difference between L1 norm of wavelets, and L1 norm given by proximal: %e.\n',L1norm-wavL1norm(wavSoftTreshold(W*x,lambdas),lambdas))
fprintf('\nIf one of the difference result was not 0, there is a problem with the implementations.\n')
