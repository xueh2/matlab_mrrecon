
%% A: Plot the k-space pattern.
close all, clear all
cd C:\MATLAB7\work\Work_Space_Dual_CPU\2011\2011_03_28_Test_3D_TGRAPPA

[Data, asc, prot] = Read_RawData('F:\Yu_Ding_Raw_Data\2011_03_16_3D_TPat\meas_MID453_CV_Perfusion_3D_TPAT_FID71270.dat');
N_1 = length(asc);
slc = zeros(14, N_1);
ssd = zeros(14, N_1);
for i=1:N_1, 
    slc(:,i)=asc(i).sLC; 
    ssd(:,i)=asc(i).sSD; 
end
figure(1), imagesc(slc)
figure(2), imagesc(ssd)

figure(1), plot(slc(1, :), 'o'), title('Phase Encoding Direction')
figure(4), plot(slc(4, :), 'o'), title('Partition Direction')
figure(8), plot(slc(7, :), 'o')

% find the first frame, and check the k-space
F_0 = find(squeeze(slc(7, :))==0);
figure(3), plot( slc(1, F_0), slc(4, F_0), 'o', 'markersize', 4, 'linewidth', 1.5)
title('{\it k}-space Acquistion Pattern', 'fontsize', 20)
xlabel('Phase Encoding Direction \rightarrow', 'fontsize', 20)
ylabel('Patition Direction \rightarrow', 'fontsize', 20, 'interpreter', 'Tex')
set(gca, 'xlim', [-1 82], 'ylim', [-0.5 11.5], 'linewidth', 2.0, 'fontsize', 16)
%print -f3 -dpdf -r600 3D_Perfusion_kspace_Acquisition_Pattern_2010_03_16 

%% B: Read in the Raw Data
close all, clear all
cd F:\Yu_Ding_Raw_Data\2011_03_16_3D_TPat
%C:\MATLAB7\work\Work_Space_Dual_CPU\2011\2011_03_28_Test_3D_TGRAPPA

[Data, asc, prot] = Read_RawData('meas_MID453_CV_Perfusion_3D_TPAT_FID71270.dat');
N_1 = length(asc);
slc = zeros(14, N_1);
for i=1:N_1, slc(:, i)=asc(i).sLC; end
%figure(1), imagesc(slc) 
Raw_Data = zeros(256, 81, 12, 12, 12, 'single');
for index = 13:length(asc)-11
    Raw_Data(:, slc(1, index)+1, slc(4, index)+1, asc(index).ushChannelId+1, slc(7, index)+1) = Data(:, index);
    % Raw_Data(fe, pe, pt, ch, temp_frames);
end
Noise_Scan = Data(:,1:12);
N_1 = length(asc)-24;
PE_Index = zeros(1, N_1);
PT_Index = zeros(1, N_1);
FR_Index = zeros(1, N_1);
j=0;
for i=13:12+N_1, 
    temp = asc(i).sLC;
    j = j + 1;
    PE_Index(:, j) = temp(1) ; 
    PT_Index(:, j) = temp(4) ; 
    FR_Index(:, j) = temp(7) ;
end
save Raw_Data_2011_03_16_3D_TPat Raw_Data Noise_Scan *_Index

%% C: Start to work on the two othorgonal 2-D kernel.
close all, clear all
load F:\Yu_Ding_Raw_Data\2011_03_16_3D_TPat\Raw_Data_2011_03_16_3D_TPat
cd C:\MATLAB7\work\Work_Space_Dual_CPU\2011\2011_03_28_Test_3D_TGRAPPA
%Raw_Data(:,:,20,:,:) = 0;
% Find the acquired pattern
F_0 = find(FR_Index==0); % First Frame
F_1 = find(FR_Index==1);
F_2 = find(FR_Index==2);
F_3 = find(FR_Index==3);
F_4 = find(FR_Index==4);
F_5 = find(FR_Index==5);
Acq_Pattern = zeros(81, 12, 12);
%Acq_Pattern_0 = zeros(81, 12);
%Acq_Pattern_1 = zeros(81, 12);
%Acq_Pattern_2 = zeros(81, 12);
%Acq_Pattern_3 = zeros(81, 12);
%Acq_Pattern_4 = zeros(81, 12);
%Acq_Pattern_5 = zeros(81, 12);
% for i=1:length(F_0)
%     Acq_Pattern_0( PE_Index(F_0(i))+1, PT_Index(F_0(i))+1 ) = 1; 
%     Acq_Pattern_1( PE_Index(F_1(i))+1, PT_Index(F_1(i))+1 ) = 1;
%     Acq_Pattern_2( PE_Index(F_2(i))+1, PT_Index(F_2(i))+1 ) = 1;
%     Acq_Pattern_3( PE_Index(F_3(i))+1, PT_Index(F_3(i))+1 ) = 1;
%     Acq_Pattern_4( PE_Index(F_4(i))+1, PT_Index(F_4(i))+1 ) = 1;
%     Acq_Pattern_5( PE_Index(F_5(i))+1, PT_Index(F_5(i))+1 ) = 1;
% end
for i=0:11
    Acq_Pattern( PE_Index( find(FR_Index==i) )+1, PT_Index( find(FR_Index==i) )+1, i+1 ) = 1; 
end
figure(3), 
% subplot('position', [0.15, 0.85, 0.81, 0.12]), imagesc(Acq_Pattern_0'), axis image, colormap(gray), axis off
% subplot('position', [0.15, 0.70, 0.81, 0.12]), imagesc(Acq_Pattern_1'), axis image, colormap(gray), axis off
% subplot('position', [0.15, 0.55, 0.81, 0.12]), imagesc(Acq_Pattern_2'), axis image, colormap(gray), axis off
% subplot('position', [0.15, 0.40, 0.81, 0.12]), imagesc(Acq_Pattern_3'), axis image, colormap(gray), axis off
% subplot('position', [0.15, 0.25, 0.81, 0.12]), imagesc(Acq_Pattern_4'), axis image, colormap(gray), axis off
% subplot('position', [0.15, 0.10, 0.81, 0.12]), imagesc(Acq_Pattern_5'), axis image, colormap(gray) 
subplot('position', [0.15, 0.85, 0.81, 0.12]), imagesc(Acq_Pattern(:,:,1)'), axis image, colormap(gray), axis off
subplot('position', [0.15, 0.70, 0.81, 0.12]), imagesc(Acq_Pattern(:,:,2)'), axis image, colormap(gray), axis off
subplot('position', [0.15, 0.55, 0.81, 0.12]), imagesc(Acq_Pattern(:,:,3)'), axis image, colormap(gray), axis off
subplot('position', [0.15, 0.40, 0.81, 0.12]), imagesc(Acq_Pattern(:,:,4)'), axis image, colormap(gray), axis off
subplot('position', [0.15, 0.25, 0.81, 0.12]), imagesc(Acq_Pattern(:,:,5)'), axis image, colormap(gray), axis off
subplot('position', [0.15, 0.10, 0.81, 0.12]), imagesc(Acq_Pattern(:,:,6)'), axis image, colormap(gray) 
text(-6, -15, 'Partition Direction', 'fontsize', 16, 'rotation', 90)
text(20, 19, 'Phase Encoding Direction', 'fontsize', 16)
%text(85, 6, '6', 'fontsize')
print -f3 -dpdf -r600 2011_03_16_3D_TPat_Acq_Pattern

figure(1), plot(PE_Index(F_0(:))+1, 'o')
figure(2), plot(PT_Index(F_0(:))+1, 'o')

%Temp = mean(Raw_Data( 64:192, 20:60, 4:8, :, : ), 5); 
Temp = mean(Raw_Data( 64:192, :, :, :, : ), 5); 
% Make sure there are no holes in the training data
sT = size(Temp); figure(4), imagesc(log(abs(reshape(Temp, sT(1), prod(sT(2:4)) ))))
Train = zeros(sT(1), sT(2), sT(3), sT(4));
Train = Temp;
%for i=1:10, Train_PE( 320*(i-1)+1:320*(i-0), :, : ) = Temp(:, :, i, :); end

acc_pe = 3;
acc_pt = 2;
temp_acq = zeros( (sT(1)-0)*(sT(2)-0)*(sT(3)), 4*4*5*sT(4) );
temp_mis = zeros( (sT(1)-0)*(sT(2)-0)*(sT(3)), (acc_pe*acc_pt-0)*sT(4)); % recon the acquired line, too
tic
shift_index = 0;
for index_pe = 1:sT(2) % Rate 3 in pe direction 
    for index_pt = 1:sT(3) % rate 2 in pt direction 
        block_pe = mod( (index_pe-3:3:(index_pe+2*3))-1, sT(2)) + 1 ; 
        block_pt = mod( (index_pt-2:2:(index_pt+2*2))-1, sT(3)) + 1 ; 
        temp = 0;
        for index_fe = 1:sT(1)
            shift_index = shift_index + 1;
            block_fe = mod( (index_fe-2:index_fe+2)-1, sT(1)) + 1;
            temp = reshape(Train(block_fe, block_pe, block_pt,:), 1, 4*4*5*12 );
            temp_acq(shift_index, :) = reshape( Train( block_fe, block_pe, block_pt, :), 1, 4*4*5*12 ); 
            mis_pe = mod( (index_pe:index_pe+(acc_pe-1))-1, sT(2) ) + 1; % missing data block
            mis_pt = mod( (index_pt:index_pt+(acc_pt-1))-1, sT(3) ) + 1; % missing data block
            temp_mis(shift_index, :) = reshape( Train(index_fe, mis_pe, mis_pt, :), 1, (acc_pe*acc_pt-0)*sT(4) ) ;   
                    
            %temp_index = ((index_pe-1)*4 + index_pt - 1)*240 + shift_index;
            %temp_acq(temp_index, :) = reshape( Train( block_fe, block_pe, block_pt, :), 1, 4*4*5*12 );
        end
        %temp_index = ((index_pe-1)*4 + index_pt -1 )*240+1:((index_pe-1)*4 + index_pt )*240 ; 
        %temp_mis(temp_index, :) = reshape( Train( 41:280, index_pe:index_pe+2, index_pt:index_pt+1, :), 240, 3*2*12 );
        %temp_mis(temp_index, :) =  ;
    end
end
toc
s_temp2 = size(temp_acq);

%%
%data_0 = (temp2' * transpose(temp1)); % I want it to be double Ding 2010-10-31
%temp2 = single(temp2); % Do not bring this to front Ding 2010-10-31
%c_0 = (temp2'*temp2); %toc, tic,'2',
%psedu_inv_1 = inv_reg(c_0, 0.0001);
C_0 = temp_acq'*temp_acq ; %+ 20*max(s_temp2(:)*var(Noise_Scan(:)))*diag(ones(1, s_temp2(2)));
dTraceR = trace(real(C_0));
nRow = size(C_0, 1);
lamda  = 0.001 * dTraceR / nRow;
%E_0 = eig( C_0 ); 
%coefficient = inv( C_0  ) * ( temp_acq' * (temp_mis) ); 
%coefficient = inv( C_0 + E_0(end)*0.0001*diag(ones(1, size(C_0, 1)))  ) * ( temp_acq' * (temp_mis) ); 
coefficient = inv( C_0 + lamda*diag(ones(1, size(C_0, 1)))  ) * ( temp_acq' * (temp_mis) ); 
tic
s_Raw = size(Raw_Data);
temp_acq_0 = zeros( prod(s_Raw(1:3))/acc_pe/acc_pt, 4*4*5*s_Raw(4), 'single' );
Temp_Recon = zeros( prod(s_Raw(1:3))/acc_pe/acc_pt, size(coefficient, 2), s_Raw(5), 'single' );
% find the first block starting point
for i=1:s_Raw(5)
    index_pe_0(i) = mod(i-1, acc_pe)+1;
    if mod(i-1, acc_pe*acc_pt)+1 < acc_pe + 1,
        index_pt_0(i) = 1;
    else
        index_pt_0(i) = 2; 
    end
end

for index_frame = 1:s_Raw(5)
    shift_index = 0;
    for index_pe = index_pe_0(index_frame):3:s_Raw(2)
        %t_i_pe = t_i_pe + 1; t_i_pt = 0;
        block_pe = mod( ((index_pe-3):3:(index_pe+3*2))-1, s_Raw(2)) + 1 ;
        for index_pt = index_pt_0(index_frame):2:s_Raw(3)
            block_pt = mod( ((index_pt-2):2:(index_pt+2*2))-1, s_Raw(3)) + 1 ;
            %shift_index = 0;
            for index_fe = 1:s_Raw(1)
                shift_index = shift_index + 1;
                block_fe = mod( (index_fe-2:index_fe+2)-1, s_Raw(1) )+1;
                temp_acq_0(shift_index, :) = reshape( Raw_Data( block_fe, block_pe, block_pt, :, index_frame), 1, 4*4*5*12 );
            end
        end
        %{index_pe, toc}
    end
    %toc
    %figure(4), imagesc(abs(temp_acq).^0.1)
    Temp_Recon(:,:,index_frame) = temp_acq_0*coefficient;
end
Recon = zeros( size(Raw_Data( :, :, :, :, :)) );
toc,
%%
% Fill k-space lines
% t_i_pe = 0; 
% for index_pe = 1:3:s_Raw(2) %2:3:80
%     t_i_pe = t_i_pe + 1; t_i_pt = 0;
%     temp = 0;
%     for index_pt = 1:2:s_Raw(3)
%         t_i_pt = t_i_pt + 1; shift_index = shift_index + 1;
%         block_pe = mod( ((index_pe-3):3:(index_pe+3*2))-1, s_Raw(2) ) +1; 
%         block_pt = mod( ((index_pt-2):2:(index_pt+2*2))-1, s_Raw(3) ) +1; 
%         %{(t_i_pe-1)*6 + t_i_pt - 1, shift_index}
%         temp_index = ((t_i_pe-1)*6 + t_i_pt - 1)*s_Raw(1) + 1 : ((t_i_pe-1)*6 + t_i_pt )*s_Raw(1);
%         temp = Temp_Recon(temp_index, :);
%         y_index = mod( (index_pe:index_pe+2)-1, s_Raw(2) )+1;
%         z_index = mod( (index_pt:index_pt+1)-1, s_Raw(3) )+1;
%         %Recon(:, index_pe:index_pe+2, index_pt:index_pt+1, :) = reshape( temp, 320, 3, 2, 12 );
%         Recon(:, y_index, z_index, :) = reshape( temp, s_Raw(1), acc_pe, acc_pt, s_Raw(4) );
%     end
%     %{index_pe, toc}
% end
for index_frame = 1:s_Raw(5)
shift_index = 0;
for index_pe = index_pe_0(index_frame):3:s_Raw(2)
    y_index = mod( (index_pe:index_pe+2)-1, s_Raw(2) )+1;
    for index_pt = index_pt_0(index_frame):2:s_Raw(3)
        z_index = mod( (index_pt:index_pt+1)-1, s_Raw(3) )+1;
        for index_fe = 1:s_Raw(1) 
            shift_index = shift_index + 1;
            Recon(index_fe, y_index, z_index, :, index_frame) = reshape( Temp_Recon(shift_index, :, index_frame), 1, acc_pe, acc_pt, s_Raw(4) );
        end
    end
    %{index_pe, toc}
end 

for i = 1:s_Raw(2)
    for j=1:s_Raw(3)
        if Acq_Pattern(i, j, index_frame) == 1,
            Recon(:,i,j,:, index_frame) = Raw_Data(:,i,j,:,index_frame);
        end
    end
end

for i=1:12,
    Recon_Image(:,:,:,i,index_frame) = ifftn( Recon(:,:,:,i, index_frame) );
    Wrap_Image(:,:,:,i, index_frame) = ifftn( Raw_Data(:,:,:,i,index_frame) );
end
end
Recon_Image = fftshift( fftshift( fftshift(Recon_Image, 3), 2), 1);
Wrap_Image = fftshift( fftshift( fftshift(Wrap_Image, 3), 2), 1);
Recon_SoS = squeeze((sqrt( mean( abs(Recon_Image).^2, 4 ) )));
Wrap_SoS = (sqrt( mean( abs(Wrap_Image).^2, 4 ) ));
close all, 
figure(1), imagesc(squeeze(abs(Recon(129,:,:,1, index_frame))).^0.1)
K_Recon_Temp = squeeze(abs(Recon(:,:,:,1, index_frame))).^0.5;
Raw_Data_Temp = squeeze(abs(Raw_Data(:,:,:,1, index_frame))).^0.5;
for i=1:12,
figure(2), subplot(2,2,1), imagesc((Recon_SoS(:,:,i,2))), title(num2str(i)), 
figure(2), subplot(2,2,2), imagesc((Wrap_SoS(:,:,i))), title(num2str(i)),
figure(2), subplot(2,2,3), imagesc(K_Recon_Temp(:,:,i)), title(num2str(i)),
figure(2), subplot(2,2,4), imagesc(Raw_Data_Temp(:,:,i)), title(num2str(i)),
pause(0.1), end

Final_Recon = Recon_SoS(s_Raw(1)/4+1:3*s_Raw(1)/4, :, :, :);
Imgs_Temp = zeros(s_Raw(1), 81*6, 1, 12 );
Imgs_Temp( 1:s_Raw(1)/2, 0*81+1:1*81, 1, 1:12 ) = Final_Recon(:,:,1,:);
Imgs_Temp( 1:s_Raw(1)/2, 1*81+1:2*81, 1, 1:12 ) = Final_Recon(:,:,2,:);
Imgs_Temp( 1:s_Raw(1)/2, 2*81+1:3*81, 1, 1:12 ) = Final_Recon(:,:,3,:);
Imgs_Temp( 1:s_Raw(1)/2, 3*81+1:4*81, 1, 1:12 ) = Final_Recon(:,:,4,:);
Imgs_Temp( 1:s_Raw(1)/2, 4*81+1:5*81, 1, 1:12 ) = Final_Recon(:,:,5,:);
Imgs_Temp( 1:s_Raw(1)/2, 5*81+1:6*81, 1, 1:12 ) = Final_Recon(:,:,6,:);

Imgs_Temp( s_Raw(1)/2+1:s_Raw(1), 0*81+1:1*81, 1, 1:12 ) = Final_Recon(:,:,7,:);
Imgs_Temp( s_Raw(1)/2+1:s_Raw(1), 1*81+1:2*81, 1, 1:12 ) = Final_Recon(:,:,8,:);
Imgs_Temp( s_Raw(1)/2+1:s_Raw(1), 2*81+1:3*81, 1, 1:12 ) = Final_Recon(:,:,9,:);
Imgs_Temp( s_Raw(1)/2+1:s_Raw(1), 3*81+1:4*81, 1, 1:12 ) = Final_Recon(:,:,10,:);
Imgs_Temp( s_Raw(1)/2+1:s_Raw(1), 4*81+1:5*81, 1, 1:12 ) = Final_Recon(:,:,11,:);
Imgs_Temp( s_Raw(1)/2+1:s_Raw(1), 5*81+1:6*81, 1, 1:12 ) = Final_Recon(:,:,12,:);

Imgs_Temp_scale = uint8(Imgs_Temp*10^5*400);
imwrite(Imgs_Temp_scale, '2011_03_16_3D_Perf_Recon_AllFrames.gif', 'DelayTime', 0.25, 'LoopCount', Inf)

figure(1), imagesc(abs(Imgs_Temp(:,:,1,6)).^0.5), axis image
% Imgs_Temp = 255*ones(s_Raw(1), s_Raw(2)*2+10, 1, 20);
% Imgs_Temp(:,  1:s_Raw(2)*1, 1, :) = reshape(Recon_SoS, s_Raw(1), s_Raw(2), 1, 20 );
% Imgs_Temp(:, 92:s_Raw(2)*2+10, 1, :) = reshape(Wrap_SoS, s_Raw(1), s_Raw(2), 1, 20 );
% Kspace_Temp = 255*ones(s_Raw(1), s_Raw(2)*2+10, 1, 20);
% Kspace_Temp(:,  1:s_Raw(2)*1, 1, :) = reshape( K_Recon_Temp, s_Raw(1), s_Raw(2), 1, 20 );
% Kspace_Temp(:, 92:s_Raw(2)*2+10, 1, :) = reshape(Raw_Data_Temp, s_Raw(1), s_Raw(2), 1, 20 );
% Imgs_Temp_uint8 = uint8(Imgs_Temp/median(Imgs_Temp(:))*6);
% %imwrite(Imgs_Temp_uint8, '2011_03_16_3D_Perf_Recon_Frame_1.gif', 'DelayTime', 0.25, 'LoopCount', Inf)
% Kspace_Temp_uint8 = uint8(Kspace_Temp/median(K_Recon_Temp(:))*64);
% %imwrite(Kspace_Temp_uint8, '2011_03_16_3D_Perf_Kspace_Frame_1.gif', 'DelayTime', 0.25, 'LoopCount', Inf)
% 


%%
x=randn(86400,960); y=randn(960,72); tic, z = x*y; toc
x=randn(86400/27,960); y=randn(960,72); tic, for i=1:27, z = x*y; end, toc

















