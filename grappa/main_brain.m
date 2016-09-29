%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Main function 
%%%%%%
%%%%%% An example for nonlinear GRAPPA method
%%%%%%
%%%%%% Written by: Yuchou Chang, University of Wisconsin - Milwaukee
%%%%%% Email: yuchou@uwm.edu; leiying@uwm.edu
%%%%%% Citation: Y. Chang, D. Liang, L. Ying, "Nonlinear GRAPPA: A Kernal Approach 
%%%%%%           to Parallel MRI Reconstruction". Magn. Reson. Med. 2012
%%%%%% Created on Oct. 12, 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
warning off all

coil_num=8;
rfac=5;%outer reduction factor

temp = zeros([256,256,coil_num]);

load rawdata_brain;

for k = 1:coil_num
    temp(:,:,k) = fftshift(ifft2(fftshift(raw_data(:,:,k))));
end

%getting the reference image
c=0;
for k=1:coil_num
    a=real(temp(:,:,k));
    b=imag(temp(:,:,k));
    c=a.*a+b.*b+c;
end
Ref_image=sqrt(c);
amp=abs(Ref_image);
amp_ref=amp;
figure;imshow(Ref_image,[0 max(amp_ref(:))*0.5]);
Img_NMSE=Ref_image/mean(mean(abs(Ref_image)));

[d1,d2,d3]=size(temp); 
ndim=d1;%phase encoding direction
off = 0;%starting sampling location

% The number of ACS lines
nencode=42;
 
% The convlution size 
num_block=2;
num_column=15;

% Obtain ACS data and undersampled data
acs_line_loc=(ndim/2+1-nencode/2):(ndim/2+nencode/2);
for l=1:coil_num
    k_space_full=fftshift(fft2(temp(:,:,l)));
    
    k_space_red(:,:,l)=k_space_full((off+1):rfac:(d1-off),:);
    acs_data(:,:,l)=k_space_full(acs_line_loc,:);
end

% Obtain uniformly undersampled locations
pe_loc=(off+1):rfac:(d1-off);

% Net reduction factor
acq_idx = zeros(d1,1);
acq_idx(pe_loc) = 1;
acq_idx(acs_line_loc) = 1;
NetR = d1 / sum(acq_idx)

% GRAPPA Reconstruction
% tic
% [full_fourier_data0, ImgRecon0, coef0] = grappa(k_space_red, rfac, pe_loc, acs_data, acs_line_loc, num_block, num_column);
% toc
% figure, imshow(abs(fftshift(ImgRecon0)),[0 max(amp_ref(:)*0.5)]);


% Nonlinear GRAPPA Reconstruction
times_comp = 3;%The number of times of the first-order terms
tic
[full_fourier_data1, ImgRecon1, coef1] = nonlinear_grappa(k_space_red, rfac, pe_loc, acs_data, acs_line_loc, num_block, num_column, times_comp);
toc
figure, imshow(abs(fftshift(ImgRecon1)),[0 max(amp_ref(:)*0.5)]);