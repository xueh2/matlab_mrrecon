
% Test 3-D SPIRiT Kernel Estimate
close all, clear all
cd C:\MATLAB7\work\Work_Space_Dual_CPU\2012\2012_09_28_Dr_Wang
load E:\Yu_Ding_Raw_Data\2012_09_28_Dr_Wang\3D_MPRAGE_Raw_Data_256_224_160_4_20120928
%Noise = Noise_Scan;
s_0 = size( K_space );

% Noise pre-whiten:
% tic
% Psi = Noise'*Noise/(size(Noise, 1)-1);
% [V, D] = eig( (Psi+Psi')/2 ); % Eigen-decomposition
% E = diag(1./sqrt(diag(D))); 
% Psi_inv_sqrt = V*E*V'; % Matrix: Psi^(-1/2) 
% temp = reshape( K_space, prod(s_0(1:3)), s_0(4) ); 
% temp = temp*Psi_inv_sqrt; % whittening 
% K_space = reshape( temp, s_0 ); 
% toc, 

figure(1), imagesc(abs(squeeze( K_space(256, :,:,1) )).^0.1)
Ref_Data = K_space( :, 102:126, :, : ); % 2: 102-126; 
s_1 = size( Ref_Data );

N_k_1 = 2; % kernel size 1 
N_k_2 = 2; % kernel size 2 
N_k_3 = 2; % kernel size 3 

N_e = (s_1(1)-2*N_k_1)*(s_1(2)-2*N_k_2)*(s_1(3)-2*N_k_3); % Number of Equations
Data_in = zeros( prod([ 2*N_k_1+1, 2*N_k_2+1, 2*N_k_3+1, s_0(4) ]) - s_0(4), N_e, 'single' ); 
Data_out = zeros( s_0(4), N_e, 'single'  ); 
Data_in = complex(Data_in, Data_in); 
Data_out = complex(Data_out, Data_out); 
% identify the sliding window pattern
temp = ones( 2*N_k_1+1, 2*N_k_2+1, 2*N_k_3+1, s_0(4) );
temp( N_k_1+1,  N_k_2+1, N_k_3+1, :) = 0;
[I, J] = find( temp(:) ==1 );

counter = 0; 
tic 
for index_1 = 1+N_k_1:s_1(1)-N_k_1
    for index_2 = 1+N_k_2:s_1(2)-N_k_2
        for index_3 = 1+N_k_3:s_1(3)-N_k_3
            counter  = counter + 1;
            temp = Ref_Data( index_1-N_k_1:index_1+N_k_1 , index_2-N_k_2:index_2+N_k_2, index_3-N_k_3:index_3+N_k_3, : ) ;
            Data_in(:, counter) = temp(I);
            Data_out(:, counter) = ( Ref_Data( index_1 , index_2, index_3, : ) );
            %Data_out(counter, :) = squeeze( K_space( index_1 , index_2, index_3, : ) );
        end
    end
end
toc  

s_D = size( Data_in ); 
tic 
M_temp = Data_in*Data_in'; 
Kernel = (Data_out*Data_in')*pinv( M_temp ); 
toc 

E_0 = eig(M_temp);
figure(2), semilogy(E_0, 'o')

disp('the std(b(:)), Ax = b')
disp(std(Data_out(:))) 
temp = Kernel*Data_in - Data_out; 
disp('std( Ax-b )')
disp(std(temp(:))) 

Diff_D = abs(Data_out(:)./temp(:));
{'mean = ', mean(Diff_D), 'median =', median(Diff_D)}
figure(3), hist(Diff_D, 1024)










