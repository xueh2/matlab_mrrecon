

function kernel = PRUNO_Kernel_2D(Ref_Img, option)
% function Coef = GRAPPA_Kernel_2D(Ref_Img, option)
% 

s_0 = size(Ref_Img);
Coef = 0;

% make sure Kernel Size Consistent
KernelSize = option.KernelSize;
Pattern_pe = option.KernelPattern;
OutPattern = option.OutPattern;

Start_fe = (KernelSize(2)-1)/2 ;
End_fe = (KernelSize(2)-1)/2 ;

Pattern_fe = -Start_fe:End_fe;

Start_pe = -min(Pattern_pe) ;
End_pe = max(Pattern_pe) ;


CoefSize = prod(KernelSize)*s_0(3);

MisSize = length(OutPattern)*s_0(3) ;
%SamplingPattern = option.SamplingPattern ;

if length(size(KernelSize))~= 2, 
    return, 
else
    Coef = zeros( CoefSize, MisSize );
end


Data_Acq = zeros( (s_0(1) - Start_fe - End_fe)*(s_0(2) - Start_pe - End_pe ), prod(KernelSize)*s_0(3), 'single' );
Data_Acq = complex( Data_Acq, Data_Acq );
%Data_Mis = zeros( (s_0(1) - Start_fe - End_fe)*(s_0(2) - Start_pe - End_pe ), length(OutPattern)*s_0(3), 'single' );
%Data_Mis = complex( Data_Mis, Data_Mis );

counter  = 0;
for index_fe = Start_fe+1:s_0(1)-End_fe
    for index_pe = Start_pe+1:s_0(2)-End_pe
        counter  = counter  + 1 ; 
        temp = Ref_Img( index_fe+Pattern_fe, index_pe+Pattern_pe, :) ;
        %keyboard
        Data_Acq(counter, :) = temp(:);
        %size(Data_Mis(counter, :)), size(Ref_Img( index_fe, index_pe+OutPattern, :))
        %Data_Mis(counter, :) = reshape( Ref_Img( index_fe, index_pe+OutPattern, :), 1, MisSize);
    end
end
counter  = counter; % test code
Asquare = (Data_Acq'*Data_Acq);
s_A = size(Asquare);
% tic
%[U,S,V] = svd(Asquare);
%kernel = U(:, s_A/2+1:end);
% toc, tic
[V, D] = eig(Asquare);
%kernel = V(:, s_A*0.75:-1:1);
%kernel = V(:, s_A-s_0(3)*1.0:-1:s_A/2+1);
kernel = V( :, s_A-s_0(3)*3.0:s_A-s_0(3)*1.0 );
E = diag(D);
figure(3), plot(cumsum(E(end:-1:1))/sum(E), 'o') 
%toc

%kernel2 = U(:, 1+32*8:s_A/2);
% close all
% figure(1), imagesc(abs(Data_Acq*kernel)), colorbar
% figure(2), imagesc(abs(Data_Acq*kernel2)), colorbar
% figure(3), imagesc(log(abs( (kernel - kernel2)./(kernel+kernel2) ))), colorbar
% figure(4), imagesc(log(abs(kernel))), colorbar
% figure(5), semilogy(diag(D(end:-1:1, end:-1:1)), 'o')











