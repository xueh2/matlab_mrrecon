% This script demonstare a CS reconstruction from 
% Randomly undersmapled phase encodes of 2D FSE
% of a brain image.

clear all
close all

addpath(strcat(pwd,'/utils'));

if exist('FWT2_PO') <2
	error('must have Wavelab installed and in the path');
end

load brain512


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(data); 	% image Size
DN = size(data); 	% data Size
TVWeight = 0.002; 	% Weight for TV penalty
xfmWeight = 0.005;	% Weight for Transform L1 penalty
Itnlim = 8;		% Number of iterations


%generate Fourier sampling operator
FT = p2DFT(mask, N, 1, 2);

% scale data
im_dc = FT'*(data.*mask./pdf);
data = data/max(abs(im_dc(:)));
im_dc = im_dc/max(abs(im_dc(:)));

%generate transform operator
XFM = Wavelet('Daubechies',4,4);	% Wavelet

% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

figure(100), imshow(abs(im_dc),[]);drawnow;

res = XFM*im_dc;

% do iterations
% tic
% for n=1:5
% 	res = fnlCg(res,param);
% 	im_res = XFM'*res;
% 	figure(100), imshow(abs(im_res),[]), drawnow
% end
% toc

% do NESTA recon
ind = find(mask>0);
b = param.data(ind(:));

U = @(z) (reshape(XFM*reshape(z, N), [N(1)*N(2) 1]));

imW = XFM*im_dc;
imW2 = U(im_dc);
norm(imW(:)-imW2)

Ut = @(z) (reshape(XFM'*reshape(z, N), [N(1)*N(2) 1]));

im = XFM'*res;
im2 = Ut(res(:));
norm(im(:)-im2)

mu = 0.01; %--- can be chosen to be small
opts = [];
opts.maxintiter = 5;
opts.TOlVar = 1e-5;
opts.verbose = 0;
opts.maxiter = 5000;
opts.U = U;
opts.Ut = Ut;
opts.stoptest = 1;  
opts.typemin = 'l1';
delta = 0;

Ac = @(z) (i2v(FT*reshape(z,N), ind));
Atc = @(z) (reshape(FT'*v2i(z,ind,N), [N(1)*N(2) 1] ));

p = Ac(im_dc);
norm(p-data(ind))

im2 = Atc(param.data(ind));
im3 = FT'*param.data;
norm(im2-im3(:))

tic;
[x_nesta,niter,resid,err] = NESTA(Ac,Atc,b,mu,delta,opts);
t.NESTA = toc
tvnesta = calctv(N(1),N(2),reshape(x_nesta,N));
im_nesta = reshape(x_nesta,N);
figure, imshow(abs(im_nesta),[]);

% show the results
figure, imshow(abs(cat(2,im_dc,im_res, im_nesta)),[]);
figure, imshow(abs(cat(2,im_dc(155:304,110:209), im_res(155:304,110:209), im_nesta(155:304,110:209))), [0,1],'InitialMagnification',200);
title(' zf-w/dc              l_1 Wavelet              l_1 Nesta');
