function kspaceCorrected = performEPIPhaseCorrection(kspace, phaseData, filterPhaseData, robustEstimation, mode)
% This function performs the EPI phase correction using OlineTSE type methods for a slice of dynamic data
% kspace : 5-D input kspace data array (Col Line Cha 1 Seg Rep/PHS)
% phaseData : 4-D phase reference data array (Col Cha Acq Seg) for this slice
% phaseLin : the line index of acquired phase line
% the SEG=1 is the regular lines; seg=2 is the reflected lines
% filterPhaseData : if 1, use a hanning window to filer the phase data
% robustEstimation : if 1, use robust estimation of mean phase slope and offset
% mode : 'PointByPoint', 'LinearPhase', 'LinearPhaseWithOffset'
% 'PointByPoint' : each point of imaging scan is corrected by multiplying the conjugate complex point in the corresponding phaseData
% 'LinearPhase' : estimate the correct the linear slope of phase changes in the imaging scans, using autocorrelation
% 'LinearPhaseWithOffset' : first estimate and correct the linear slope of phase changes, then correct offset of phases
% [Col Line Cha Acq Slice Partition Echo Phase Rep Set Seg]

% perform inverse fft along the fe direction
% kspace = fftshift(kspace,1);
% kspace = ifftshift(ifft(kspace, [], 1), 1);

kspace = ifftc(kspace, 1);
ss = size(kspace);

if ( length(ss)==5 )
    ss = [ss 1];
end
% 
% kspace = reshape(kspace, [ss(1) ss(2) ss(3) 1 ss(4) ss(5)]);
% ss = size(kspace);

N = ss(1);

% check whether the SEG 0 and 1 shall be swapped
swapSEG = 0;
if ( numel(size(phaseData)) == 4 )
    sumV = sum(abs(phaseData(:,1,2,2,1)));
else
    sumV = 0;
end

if ( sumV > 0 )
    swapSEG = 1;
    
    kspace2 = kspace;
    kspace2(:,:,:,1,1,:) = kspace(:,:,:,1,2,:);
    kspace2(:,:,:,1,2,:) = kspace(:,:,:,1,1,:);
    kspace = kspace2;
    
    phaseData2 = phaseData;
    phaseData2(:,:,:,1,:) = phaseData(:,:,:,2,:);
    phaseData2(:,:,:,2,:) = phaseData(:,:,:,1,:);
    phaseData = phaseData2;
end

% phaseData = fftshift(phaseData,1);
% phaseData = ifftshift(ifft(phaseData, [], 1), 1);

phaseData = ifftc(phaseData, 1);

if ( numel(size(phaseData)) == 4 ) % use ACQ and SEG to store phase corr lines
    % positive signal
    S1 = phaseData(:,:,1,1);
    S3 = phaseData(:,:,2,1);
    % compensate for different TEs
    S2 = (S1+S3)/2;
    S2 = squeeze(S2); % COL, CHA

    % negative signal
    S2R = phaseData(:,:,1,2); % reflected lines
    S2R = squeeze(S2R); % COL, CHA
else
    % positive signal
    S1 = phaseData(:,:,1);
    if ( size(phaseData, 3) >= 3 )
        S3 = phaseData(:,:,3);
        % compensate for different TEs
        S2 = (S1+S3)/2;
        S2 = squeeze(S2); % COL, CHA
    else
        S2 = squeeze(S1);
    end
    % negative signal
    S2R = phaseData(:,:,2); % reflected lines
    S2R = squeeze(S2R); % COL, CHA
end

% if required, filering the phase data
if ( filterPhaseData )
    
%     S2 = fftshift(S2,1);
%     S2 = ifftshift(fft(S2, [], 1), 1);
%     
%     S2R = fftshift(S2R,1);
%     S2R = ifftshift(fft(S2R, [], 1), 1);
    
    iWinLen = 8;
    iLowerBorder = iWinLen/2;
    iUpperBorder = N-iWinLen/2;
    twoPi  = 2. * pi / iWinLen;
    
    for cha=1:ss(3)
        for col=0:N-1            
            if ( (col<iLowerBorder) | (col>iUpperBorder) )
                S2(col+1,cha) = 0;
                S2R(col+1,cha) = 0;
            else
                filterVal = 0.5+0.5*cos(twoPi*(col-N/2));
                S2(col+1,cha) = S2(col+1,cha)*filterVal;
                S2R(col+1,cha) = S2R(col+1,cha)*filterVal;
            end            
        end
    end
    
%     S2 = fftshift(S2,1);
%     S2 = ifftshift(ifft(S2, [], 1), 1);
%     
%     S2R = fftshift(S2R,1);
%     S2R = ifftshift(ifft(S2R, [], 1), 1);    
end

% compute the linear phase changes
estimationWindow = 1;
range = [1+estimationWindow:N-estimationWindow];

if ( robustEstimation )
    
    % positive scan
    % compute the auto-correlation
    S2Sum = sum(abs(S2), 2);
    indValid = find(S2Sum>0);
    autoCorrValue = S2(range,:).*conj(S2(range-1,:));    
    autoCorrValue = autoCorrValue(indValid, :);
    % compute the angle
    angleValue = angle(autoCorrValue);
    % compute the robust mean
    options.lts = 1;
    options.alpha = 0.5;
    res.center = zeros(1,ss(3));
    for cha=1:ss(3)        
%         [resCha,raw]=fastmcd(double(angleValue(:,cha)),options);
%         res.center(cha) = resCha.center;
        phi = robustfit(range(indValid), cumsum(angleValue(:,cha)),'bisquare',4.685, 'off');
        res.center(cha) = phi(1);
    end
    
    deltaPhi = cos(res.center) + i*sin(res.center);
    deltaPhiLine = zeros(N,ss(3),ss(4)); % COL CHA SLC
    for kk=1:N
        tt = conj(deltaPhi).^(kk-1-N/2);
        for slc=1:ss(4)
            for cha=1:ss(3)
                deltaPhiLine(kk,cha,slc) = tt(cha);
            end
        end
    end
    
    % negative scan
    % compute the auto-correlation
    S2RSum = sum(abs(S2R), 2);
    indValid = find(S2RSum>0);
    autoCorrValue = S2R(range,:).*conj(S2R(range-1,:));    
    autoCorrValue = autoCorrValue(indValid, :);
    % compute the angle
    angleValue = angle(autoCorrValue);
    % compute the robust mean
%     [res,raw]=fastmcd(double(angleValue),options);
    res.center = zeros(1,ss(3));
%     for cha=1:ss(3)        
%         [resCha,raw]=fastmcd(double(angleValue(:,cha)),options);
%         res.center(cha) = resCha.center;
%     end
    for cha=1:ss(3)        
%         [resCha,raw]=fastmcd(double(angleValue(:,cha)),options);
%         res.center(cha) = resCha.center;
        phi = robustfit(range(indValid), cumsum(angleValue(:,cha)),'bisquare',4.685, 'off');
        res.center(cha) = phi(1);
    end
    
    deltaPhiR = cos(res.center) + i*sin(res.center);
    deltaPhiRLine = zeros(N,ss(3),ss(4)); % COL CHA SLC
    for kk=1:N
        tt = conj(deltaPhiR).^(kk-1-N/2);
        for slc=1:ss(4)
            for cha=1:ss(3)
                deltaPhiRLine(kk,cha,slc) = tt(cha);
            end
        end
    end
    
else
    sSum = squeeze(sum(S2(range,:).*conj(S2(range-1,:)), 1)); % 1 CHA
    deltaPhi = sSum./abs(sSum);
    deltaPhiLine = zeros(N,ss(3),ss(4)); % COL CHA SLC
    for kk=1:N
        tt = conj(deltaPhi).^(kk-1-N/2);
        for slc=1:ss(4)
            for cha=1:ss(3)
                deltaPhiLine(kk,cha,slc) = tt(cha);
            end
        end
    end

    sSum = squeeze(sum(S2R(range,:).*conj(S2R(range-1,:)), 1)); % 1 CHA
    deltaPhiR = sSum./abs(sSum);
    deltaPhiRLine = zeros(N,ss(3),ss(4)); % COL CHA SLC
    for kk=1:N
        tt = conj(deltaPhiR).^(kk-1-N/2);
        for slc=1:ss(4)
            for cha=1:ss(3)
                deltaPhiRLine(kk,cha,slc) = tt(cha);
            end
        end
    end

%% test

%     sSum = squeeze(sum(S2(range,:).*conj(S2(range-1,:)), 1)); % 1 CHA
%     angleValue = angle(sSum./abs(sSum));
%     
%     deltaPhi = cos(angleValue) + i*sin(angleValue);
%     deltaPhiLine = zeros(N,ss(3),ss(4)); % COL CHA SLC
%     for kk=1:N
%         tt = conj(deltaPhi).^(kk-1-N/2);
%         for slc=1:ss(4)
%             for cha=1:ss(3)
%                 deltaPhiLine(kk,cha,slc) = tt(cha);
%             end
%         end
%     end
% 
%     sSum = squeeze(sum(S2R(range,:).*conj(S2R(range-1,:)), 1)); % 1 CHA
%     angleValue = angle(sSum./abs(sSum));
%     
%     deltaPhiR = cos(angleValue) + i*sin(angleValue);
%     deltaPhiRLine = zeros(N,ss(3),ss(4)); % COL CHA SLC
%     for kk=1:N
%         tt = conj(deltaPhiR).^(kk-1-N/2);
%         for slc=1:ss(4)
%             for cha=1:ss(3)
%                 deltaPhiRLine(kk,cha,slc) = tt(cha);
%             end
%         end
%     end
    
%     % positive scan
%     % compute the auto-correlation
%     S2Sum = sum(abs(S2), 2);
%     indValid = find(S2Sum>0);
%     autoCorrValue = S2(range,:).*conj(S2(range-1,:));    
%     autoCorrValue = autoCorrValue(indValid, :);
%     % compute the angle
%     autoCorrValue = sum(autoCorrValue,1);
%     angleValue = angle(autoCorrValue);
%     % compute the robust mean
% %     res.center = mean(angleValue);
%     res.center = angleValue;
%     
%     deltaPhi = cos(res.center) + i*sin(res.center);
%     deltaPhiLine = zeros(N,ss(3),ss(4)); % COL CHA SLC
%     for kk=1:N
%         tt = conj(deltaPhi).^(kk-1-N/2);
%         for slc=1:ss(4)
%             for cha=1:ss(3)
%                 deltaPhiLine(kk,cha,slc) = tt(cha);
%             end
%         end
%     end
%     
%     % negative scan
%     % compute the auto-correlation
%     S2RSum = sum(abs(S2R), 2);
%     indValid = find(S2RSum>0);
%     autoCorrValue = S2R(range,:).*conj(S2R(range-1,:));    
%     autoCorrValue = autoCorrValue(indValid, :);
%     autoCorrValue = sum(autoCorrValue,1);
%     % compute the angle
%     angleValue = angle(autoCorrValue);
%     % compute the robust mean
% %     res.center = mean(angleValue);
%     res.center = angleValue;
%     
%     deltaPhiR = cos(res.center) + i*sin(res.center);
%     deltaPhiRLine = zeros(N,ss(3),ss(4)); % COL CHA SLC
%     for kk=1:N
%         tt = conj(deltaPhiR).^(kk-1-N/2);
%         for slc=1:ss(4)
%             for cha=1:ss(3)
%                 deltaPhiRLine(kk,cha,slc) = tt(cha);
%             end
%         end
%     end

end

% 'PointByPoint'
if ( strcmp(mode, 'PointByPoint') )
    
    S2Line = zeros(N,ss(3),ss(4)); % COL CHA SLC
    for slc=1:ss(4)
        S2Line(:,:,slc) = conj(S2)./(abs(S2)+eps);
    end
    
    S2RLine = zeros(N,ss(3),ss(4)); % COL CHA SLC
    for slc=1:ss(4)
        S2RLine(:,:,slc) = conj(S2R)./(abs(S2R)+eps);
    end

    for rep=1:ss(6)
        for lin=1:ss(2)
            kspace(:,lin,:,:,1,rep) = kspace(:,lin,:,:,1,rep).*reshape(S2Line,[ss(1) 1 ss(3) ss(4)]); % seg 1
            kspace(:,lin,:,:,2,rep) = kspace(:,lin,:,:,2,rep).*reshape(S2RLine,[ss(1) 1 ss(3) ss(4)]); % seg 2
        end
    end
end

% 'LinearPhase'
if ( strcmp(mode, 'LinearPhase') )
    for rep=1:ss(6)
        for lin=1:ss(2)
            kspace(:,lin,:,:,1,rep) = kspace(:,lin,:,:,1,rep).*reshape(deltaPhiLine,[ss(1) 1 ss(3) ss(4)]); % seg 1
            kspace(:,lin,:,:,2,rep) = kspace(:,lin,:,:,2,rep).*reshape(deltaPhiRLine,[ss(1) 1 ss(3) ss(4)]); % seg 2
        end
    end
end

% 'LinearPhaseWithOffset'
if ( strcmp(mode, 'LinearPhaseWithOffset') )
    
    % first correct linear phase changes
    for rep=1:ss(6)
        for lin=1:ss(2)
            kspace(:,lin,:,:,1,rep) = kspace(:,lin,:,:,1,rep).*reshape(deltaPhiLine,[ss(1) 1 ss(3) ss(4)]); % seg 1
            kspace(:,lin,:,:,2,rep) = kspace(:,lin,:,:,2,rep).*reshape(deltaPhiRLine,[ss(1) 1 ss(3) ss(4)]); % seg 2
        end
    end
    
    % compute offset phase
    S2Corrected = S2.*deltaPhiLine(:,:,1);
    S2RCorrected = S2R .*deltaPhiRLine(:,:,1);
    
    sSum = sum(S2Corrected.*conj(S2RCorrected), 1); % 1 CHA
    Phi0 = squeeze(sSum./abs(sSum));      
    Phi0Line = zeros(N,ss(3),ss(4));
    for kk=1:N
        for slc=1:ss(4)
            for cha=1:ss(3)
                Phi0Line(kk,cha,slc) = Phi0(cha);
            end
        end
    end

    for rep=1:ss(6)
        for lin=1:ss(2)
            kspace(:,lin,:,:,2,rep) = kspace(:,lin,:,:,2,rep).*reshape(Phi0Line,[ss(1) 1 ss(3) ss(4)]); % seg 2
        end
    end
end

% go back to kspace
kspace = squeeze(kspace);

if ( swapSEG )   
    kspace2 = kspace;
    kspace2(:,:,:,1,:) = kspace(:,:,:,2,:);
    kspace2(:,:,:,2,:) = kspace(:,:,:,1,:);
    kspace = kspace2;
end

% kspace = fftshift(kspace,1);
% kspaceCorrected = ifftshift(fft(kspace, [], 1), 1);
kspaceCorrected = fftc(kspace, 1);
