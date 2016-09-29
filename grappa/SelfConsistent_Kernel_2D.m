

function Coef = SelfConsistent_Kernel_2D(Ref_Img, option)
% function Coef = GRAPPA_Kernel_2D(Ref_Img, option)
% This is to use all the acquired k-space points in the sliding window as
% output, unlike GRAPPA kernel, only use the center points as output.
% Therefore, this is very close to the null-space with modelling method.
% In rate 4 data, if using two blocks, then the # of equations = 1/6 # of
% unknown, independent of the width of the 2-D kernel.

% Yu Ding 2011-07-13

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

MisSize = length(OutPattern)*KernelSize(2)*s_0(3) ;
%SamplingPattern = option.SamplingPattern ;

if length(size(KernelSize))~= 2, 
    return, 
else
    Coef = zeros( CoefSize, MisSize );
end


Data_Acq = zeros( (s_0(1) - Start_fe - End_fe)*(s_0(2) - Start_pe - End_pe ), prod(KernelSize)*s_0(3), 'single' );
Data_Acq = complex( Data_Acq, Data_Acq );
%Data_Mis = zeros( (s_0(1) - Start_fe - End_fe)*(s_0(2) - Start_pe - End_pe ), length(OutPattern)*s_0(3), 'single' );
Data_Mis = zeros( (s_0(1) - Start_fe - End_fe)*(s_0(2) - Start_pe - End_pe ), MisSize, 'single' );
Data_Mis = complex( Data_Mis, Data_Mis );

counter  = 0;
for index_fe = Start_fe+1:s_0(1)-End_fe
    for index_pe = Start_pe+1:s_0(2)-End_pe
        counter  = counter  + 1 ; 
        temp = Ref_Img( index_fe+Pattern_fe, index_pe+Pattern_pe, :) ;
        %keyboard
        Data_Acq(counter, :) = temp(:);
        %size(Data_Mis(counter, :)), size(Ref_Img( index_fe+Pattern_fe, index_pe+OutPattern, :))
        Data_Mis(counter, :) = reshape( Ref_Img( index_fe+Pattern_fe, index_pe+OutPattern, :), 1, MisSize);
    end
end
counter  = counter;
Asquare = (Data_Acq'*Data_Acq);
Tr_A = trace(Asquare) ;
%0.0001*Tr_A/CoefSize;
Coef = inv(Asquare + 0.001*Tr_A/CoefSize*eye(CoefSize) ) * (Data_Acq' * Data_Mis ) ;


































