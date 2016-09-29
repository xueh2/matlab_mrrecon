% this is a script to demonstrate the original experiment by Candes, Romberg and Tao
%
% (c) Michael Lustig 2007

rand('twister',2000);
addpath(strcat(pwd,'/utils'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = [256,256]; 		% image Size
DN = [256,256]; 	% data Size
pctg = [0.33];  	% undersampling factor
P = 5;			% Variable density polymonial degree
TVWeight = 0.01; 	% Weight for TV penalty
xfmWeight = 0.00;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations


% generate variable density random sampling
pdf = genPDF(DN,P,pctg , 2 ,0.1,0);	% generates the sampling PDF
k = genSampling(pdf,10,60);		% generates a sampling pattern

%generate image
im = (phantom(N(1)))  + randn(N)*0.01 + i*randn(N)*0.01;

%generate Fourier sampling operator
FT = p2DFT(k, N, 1, 2);
data = FT*im;

%generate transform operator

%XFM = Wavelet('Daubechies',6,4);	% Wavelet
%XFM = TIDCT(8,4);			% DCT
XFM = 1;				% Identity transform 	

% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

im_dc = FT'*(data./pdf);	% init with zf-w/dc (zero-fill with density compensation)
figure(100), imshow(abs(im_dc),[]);drawnow;

res = XFM*im_dc;

% do iterations
tic
for n=1:8
	res = fnlCg(res,param);
	im_res = XFM'*res;
	figure(100), imshow(abs(im_res),[]), drawnow
end
toc

% run the NESTA optimization
ind = find(k>0);
b = param.data(ind(:));

U = @(z) z;
Ut = @(z) z;
mu = 0.01; %--- can be chosen to be small
opts = [];
opts.maxintiter = 5;
opts.TOlVar = 1e-5;
opts.verbose = 0;
opts.maxiter = 5000;
opts.U = U;
opts.Ut = Ut;
opts.stoptest = 1;  
opts.typemin = 'tv';
delta = 0;

Ac = @(z) (i2v(FT*reshape(z,N), ind));
Atc = @(z) (reshape(FT'*v2i(z,ind,N), [N(1)*N(2) 1] ));

p = Ac(im);
norm(p-data(ind))

im2 = Atc(data(ind));
im3 = FT'*data;
norm(im2-im3(:))

tic;
[x_nesta,niter,resid,err] = NESTA(Ac,Atc,b,mu,delta,opts);
t.NESTA = toc
NA_nesta = counter();
NA_nesta_chg = 2*niter; %-- with the change of variable simplification, NA = 2*niter
tvnesta = calctv(N(1),N(2),reshape(x_nesta,N));
Xnesta = reshape(x_nesta,N);
figure, imshow(abs(Xnesta),[]);

% create a low-res mask
mask_lr = genLRSampling_pctg(DN,pctg,1,0);
im_lr = ifft2c(zpad(fft2c(im).*mask_lr,N(1),N(2)));

im_full = ifft2c(zpad(fft2c(im),N(1),N(2)));
figure, imshow(abs(cat(2,im_full,im_lr,im_dc,im_res,Xnesta)),[]);
title('original             low-res              zf-w/dc              TV        Nesta');

figure, plot(1:N(1), abs(im_full(end/2,:)),1:N(1), abs(im_lr(end/2,:)), 1:N(2), abs(im_dc(end/2,:)), 1:N(2), abs(im_res(end/2,:)), 1:N(2), abs(Xnesta(end/2,:)), 'LineWidth',2);
legend('original', 'LR', 'zf-w/dc', 'TV', 'Nesta');


