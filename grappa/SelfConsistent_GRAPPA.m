
function [Kspace, Kspace0] = SelfConsistent_GRAPPA(Ref_Img, Reduced_Kspace, GRAPPA_option, SelfConsist_option)

% %tic
% Coef = GRAPPA_Kernel_2D(Ref_Img, GRAPPA_option); 
% %figure(1), imagesc(abs(Coef)), title('1. GRAPPA Coef'), pause, 
% %toc, tic, 
% SelfConsist_Coef = GRAPPA_Kernel_2D(Ref_Img, SelfConsist_option);
% %PRUNO_kernel = PRUNO_Kernel_2D(Ref_Img, SelfConsist_option);
% %figure(1), imagesc(abs(SelfConsist_Coef)), title('2. SelfConsist Coef'), pause, 
% %toc, tic, 
% Kspace0 = GRAPPA_Recon_2D(Reduced_Kspace, Coef, GRAPPA_option);
% %figure(1), imagesc(abs(Kspace0(:,:,1))), title('3. GRAPPA Kspace'), pause, 
% %toc, tic, 
% Kspace = SelfConsistent_GRAPPA_Recon_2D(Kspace0, SelfConsist_Coef, SelfConsist_option);
% %toc


% Coef = GRAPPA_Kernel_2D(Ref_Img, GRAPPA_option); 
% %figure(1), imagesc(abs(Coef)), title('1. GRAPPA Coef'), pause, 
% %toc, tic, 
% PRUNO_kernel = PRUNO_Kernel_2D(Ref_Img, SelfConsist_option);
% size(PRUNO_kernel)
% %figure(1), imagesc(abs(SelfConsist_Coef)), title('2. SelfConsist Coef'), pause, 
% %toc, tic, 
% Kspace0 = GRAPPA_Recon_2D(Reduced_Kspace, Coef, GRAPPA_option);
% %figure(1), imagesc(abs(Kspace0(:,:,1))), title('3. GRAPPA Kspace'), pause, 
% %toc, tic, 
% Kspace = PRUNO_SelfConsistent_GRAPPA_Recon_2D(Kspace0, PRUNO_kernel, SelfConsist_option);
% %toc
% %figure(10), imagesc(abs([(Kspace(:,:,1)-Kspace0(:,:,1))*10, Kspace(:,:,1)]).^0.5), title('4. SelfConsistent GRAPPA Kspace'),pause(0.1), 


Coef = GRAPPA_Kernel_2D(Ref_Img, GRAPPA_option); 
%figure(1), imagesc(abs(Coef)), title('1. GRAPPA Coef'), pause, 
%toc, tic, 
SelfConsist_Coef = GRAPPA_Kernel_2D(Ref_Img, SelfConsist_option);
PRUNO_kernel = PRUNO_Kernel_2D(Ref_Img, SelfConsist_option);
disp(['PRUNO_kernel size', num2str(size(PRUNO_kernel)), 'SelfConsist_Coef size', num2str(size(SelfConsist_Coef))])
%figure(1), imagesc(abs(SelfConsist_Coef)), title('2. SelfConsist Coef'), pause, 
%toc, tic, 
Kspace0 = GRAPPA_Recon_2D(Reduced_Kspace, Coef, GRAPPA_option);
%figure(1), imagesc(abs(Kspace0(:,:,1))), title('3. GRAPPA Kspace'), pause, 
%toc, tic, 
Kspace = PRUNO_SelfConsistent_GRAPPA_Recon_2D(Kspace0, PRUNO_kernel, SelfConsist_option);
%toc
%figure(10), imagesc(abs([(Kspace(:,:,1)-Kspace0(:,:,1))*10, Kspace(:,:,1)]).^0.5), title('4. SelfConsistent GRAPPA Kspace'),pause(0.1), 





