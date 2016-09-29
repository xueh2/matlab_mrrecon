

function Coef = GRAPPA_Kernel_2D_2(kspace, acs, option, thresReg)
% function Coef = GRAPPA_Kernel_2D_2(kspace, acs, option)
% option.acqPELines, the index of acquired PE lines in kspace, starting from 1
% option.acqAcsLines, the index of acquired acs lines in the coordinate of kspace

s_0 = size(kspace);
Coef = 0;

Npe = size(kspace, 2);
NpeAcs = size(acs, 2);

acqPELines = option.acqPELines;
acqAcsLines = option.acqAcsLines;

% make sure Kernel Size Consistent
KernelSize = option.KernelSize;
Pattern_pe = option.KernelPattern;
OutPattern = option.OutPattern;

% the line offsets for the block patterns to every output lines
lineOffsets = zeros(numel(option.KernelPattern), numel(OutPattern));
for o=1:numel(OutPattern)
    lineOffsets(:,o) = option.KernelPattern - OutPattern(o);
end

Start_fe = (KernelSize(2)-1)/2 ;
End_fe = (KernelSize(2)-1)/2 ;

Pattern_fe = -Start_fe:End_fe;

Start_pe = -min(Pattern_pe) ;
End_pe = max(Pattern_pe) ;

CoefSize = prod(KernelSize)*s_0(3);
MisSize = length(OutPattern)*s_0(3) ;

if length(size(KernelSize))~= 2, 
    return, 
else
    Coef = zeros( CoefSize, MisSize );
end

Data_Acq = zeros( (s_0(1) - Start_fe - End_fe)*NpeAcs, prod(KernelSize)*s_0(3), 'single' );
Data_Acq = complex( Data_Acq, Data_Acq );

Data_Mis = zeros( (s_0(1) - Start_fe - End_fe)*NpeAcs, length(OutPattern)*s_0(3), 'single' );
Data_Mis = complex( Data_Mis, Data_Mis );

counter  = 0;
for index_pe = 1:NpeAcs        
    
    % blockNeeded = Pattern_pe+option.acqAcsLines(index_pe);
    blockNeeded = Pattern_pe+option.acqAcsLines(index_pe)-1;
    ind = find(blockNeeded<=0);
    blockNeeded(ind) = blockNeeded(ind) + Npe;
    ind = find(blockNeeded>Npe);
    blockNeeded(ind) = blockNeeded(ind) - Npe; 

    allAcquired = isempty(setdiff(blockNeeded, option.acqPELines))==1;        
    if ( ~allAcquired )
        continue;
    end
    
    acsNeeded = index_pe+OutPattern-1;
    ind = find(acsNeeded<=0);
    acsNeeded(ind) = acsNeeded(ind) + NpeAcs;
    ind = find(acsNeeded>NpeAcs);
    acsNeeded(ind) = acsNeeded(ind) - NpeAcs;

    for index_fe = Start_fe+1:s_0(1)-End_fe
                             
        counter  = counter  + 1; 

        temp = kspace( index_fe+Pattern_fe, blockNeeded, :);
        Data_Acq(counter, :) = temp(:);
                    
        Data_Mis(counter, :) = reshape( acs( index_fe, acsNeeded, :), 1, MisSize);
    end
end

Data_Acq = Data_Acq(1:counter, :);
Data_Mis = Data_Mis(1:counter, :);

Asquare = (Data_Acq'*Data_Acq);
AsquareInv = inv_tikhonov_IcePAT(Asquare, thresReg);
Coef = AsquareInv * (Data_Acq'*Data_Mis);
