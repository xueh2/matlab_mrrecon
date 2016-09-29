
function Kspace = PRUNO_SelfConsistent_GRAPPA_Recon_2D(Kspace, PRUNO_Coef, SelfConsist_option)

s_0 = size(Kspace);

KernelSize = SelfConsist_option.KernelSize;
Pattern_pe = SelfConsist_option.KernelPattern;
OutPattern = SelfConsist_option.OutPattern;
% Check the SelfConsist_option.OutPattern
~(OutPattern);
length(OutPattern)==1;
if (OutPattern)
    disp('SelfConsist_option.OutPattern ~=0, Error!'),
    return
elseif ~(length(OutPattern)==1)
    disp('length(SelfConsist_option.OutPattern) ~=1, Error!'),
    return
end

Start_fe = (KernelSize(2)-1)/2 ;
End_fe = (KernelSize(2)-1)/2 ;

Pattern_fe = -Start_fe:End_fe;

Start_pe = -min(Pattern_pe) ;
End_pe = max(Pattern_pe) ;


CoefSize = prod(KernelSize)*s_0(3);

MisSize = length(OutPattern)*s_0(3) ;
OutSize = length(OutPattern);
SamplingPattern = double(SelfConsist_option.SamplingPattern) ;

if length(size(KernelSize))~= 2, 
    return, 
end

Data_Acq = zeros( s_0(1)*length(SamplingPattern), s_0(3), 'single' );
Data_Acq = complex( Data_Acq, Data_Acq );
Data_Mis = zeros( s_0(1)*length(SamplingPattern), prod(KernelSize)*s_0(3), 'single' );
Data_Mis = complex( Data_Mis, Data_Mis );

counter  = 0;
for index_fe = 1:s_0(1)
    for index_pe = 1:length(SamplingPattern)
        counter  = counter  + 1 ;
        k_fe = mod(index_fe+Pattern_fe-1, s_0(1) ) + 1 ;
        k_pe = mod(SamplingPattern(index_pe) + Pattern_pe - 1, s_0(2) ) + 1;
        temp = Kspace( k_fe, k_pe, :) ;
        Data_Mis(counter, :) = temp(:);
        Data_Acq(counter, :) = Kspace( index_fe, SamplingPattern(index_pe) + OutPattern, :);
    end
end

% Reconstruction:  R*A'*inv(Q+A*R*A')
%R = 1.0;% 2.5;
%tic
%A = R*inv(eye(size(PRUNO_Coef, 2)) + R*PRUNO_Coef'*PRUNO_Coef)*PRUNO_Coef';
A = PRUNO_Coef(:,:,2)';
%{'1', toc}
%A = 1*inv( 0 + 1*PRUNO_Coef'*PRUNO_Coef)*PRUNO_Coef';
%'size A = ', svd(A), % test code Ding 2011-07-11
%'size A = ', size(A)
%'size (Data_Acq - Data_Mis*PRUNO_Coef) = ', size((Data_Acq - Data_Mis*PRUNO_Coef))
%'SVD (Data_Acq - Data_Mis*PRUNO_Coef) = ', svd((Data_Acq - Data_Mis*PRUNO_Coef))'
%'SVD ( Data_Mis) = ', svd( Data_Mis )'

% Reconstruction: x = x_0 + R*A'*inv(Q+A*R*A')*(b-A*x_0)
%tic
% temp = 0; 
% %temp = (Data_Acq - Data_Mis*PRUNO_Coef);
% temp = (0 - Data_Mis*PRUNO_Coef);
% temp1 = temp*A;
% {'Diff E.', std(temp(:)), 'Acq. E.', std(Data_Acq(:)), 'Mis. E.', std(Data_Mis(:)), 'Corr. E. ', std(temp1(:)) }
%Data_Mis = Data_Mis + (Data_Acq - Data_Mis*PRUNO_Coef) * A ;
%{'2', toc}
%tic
Data_Mis = Data_Mis + (0 - Data_Mis*PRUNO_Coef(:,:,1)) * A ;
%{'3', toc}

%Data_Mis = Data_Mis ; % Test code
%tic
% Reshape to fill k-space
counter = 0; temp = zeros(size(Kspace));
temp( :, SamplingPattern, :) = Kspace( :, SamplingPattern, :); 
Kspace = temp ;
temp = 0;
for index_fe = 1:s_0(1)
    for index_pe = 1:length(SamplingPattern)
        counter  = counter  + 1 ;
        k_fe = mod(index_fe+Pattern_fe-1, s_0(1) ) + 1 ;
        k_pe = mod(SamplingPattern(index_pe) + Pattern_pe - 1, s_0(2) ) + 1;
        % The SelfConsist_option.Averages is how many times the k-space is updated.
        temp = reshape( Data_Mis(counter, :), length(k_fe), length(k_pe), s_0(3) ) / SelfConsist_option.Averages ;
        Kspace( k_fe, k_pe, :) = Kspace( k_fe, k_pe, :) + temp; 
    end
end
%{'4', toc}




