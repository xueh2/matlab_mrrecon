function [kspace, Noise, ref] = Parse_ICE_Rawdata_VD11A_Jun(measData, flagUndoOS, flagIncludeSeparateRef)
% Data: k-space lines read from meas file, to be filled in the right
% k-space coordinate
% asc: mdf for each line
% prot: header from meas file in double format
% flagUndoOS: undo oversampling in the FE direction
% flagIncludeSeparateRef: include separate reference line for recon or not
%
% fix a bug with regard to totalLineCount on January 16, 2012
%
% The code assumes partition=1

% read the Ice meas data
% /*--------------------------------------------------------------------------*/
% /* Definition of loop counter structure                                     */
% /* Note: any changes of this structure affect the corresponding swapping    */
% /*       method of the measurement data header proxy class (MdhProxy)       */
% /*--------------------------------------------------------------------------*/
% typedef struct
% {
%   PACKED_MEMBER( uint16_t,  ushLine         ); /* line index                   */
%   PACKED_MEMBER( uint16_t,  ushAcquisition  ); /* acquisition index            */
%   PACKED_MEMBER( uint16_t,  ushSlice        ); /* slice index                  */
%   PACKED_MEMBER( uint16_t,  ushPartition    ); /* partition index              */
%   PACKED_MEMBER( uint16_t,  ushEcho         ); /* echo index                   */
%   PACKED_MEMBER( uint16_t,  ushPhase        ); /* phase index                  */
%   PACKED_MEMBER( uint16_t,  ushRepetition   ); /* measurement repeat index     */
%   PACKED_MEMBER( uint16_t,  ushSet          ); /* set index                    */
%   PACKED_MEMBER( uint16_t,  ushSeg          ); /* segment index  (for TSE)     */
%   PACKED_MEMBER( uint16_t,  ushIda          ); /* IceDimension a index         */
%   PACKED_MEMBER( uint16_t,  ushIdb          ); /* IceDimension b index         */
%   PACKED_MEMBER( uint16_t,  ushIdc          ); /* IceDimension c index         */
%   PACKED_MEMBER( uint16_t,  ushIdd          ); /* IceDimension d index         */
%   PACKED_MEMBER( uint16_t,  ushIde          ); /* IceDimension e index         */
% } sLoopCounter;
% /*--------------------------------------------------------------------------*/
% /*  Definition of measurement data header                                   */
% /*--------------------------------------------------------------------------*/
% typedef struct
% {
%   PACKED_MEMBER( uint32_t,     ulFlagsAndDMALength           );    // bit  0..24: DMA length [bytes]
%                                                                    // bit     25: pack bit
%                                                                    // bit 26..31: pci_rx enable flags                   4 byte
%   PACKED_MEMBER( int32_t,      lMeasUID                      );    // measurement user ID                               4
%   PACKED_MEMBER( uint32_t,     ulScanCounter                 );    // scan counter [1...]                               4
%   PACKED_MEMBER( uint32_t,     ulTimeStamp                   );    // time stamp [2.5 ms ticks since 00:00]             4
%   PACKED_MEMBER( uint32_t,     ulPMUTimeStamp                );    // PMU time stamp [2.5 ms ticks since last trigger]  4
%   PACKED_MEMBER( uint32_t,     aulEvalInfoMask[MDH_NUMBEROFEVALINFOMASK] ); // evaluation info mask field           8
%   PACKED_MEMBER( uint16_t,     ushSamplesInScan              );    // # of samples acquired in scan                     2
%   PACKED_MEMBER( uint16_t,     ushUsedChannels               );    // # of channels used in scan                        2   =32
%   sLoopCounter sLC;                                                                                                // loop counters                                    28   =60
%   sCutOffData sCutOff;                                                                                     // cut-off values                                    4
%   PACKED_MEMBER( uint16_t,     ushKSpaceCentreColumn         );    // centre of echo                                    2
%   PACKED_MEMBER( uint16_t,     ushCoilSelect                 );    // Bit 0..3: CoilSelect                              2
%   PACKED_MEMBER( float,        fReadOutOffcentre             );    // ReadOut offcenter value                           4
%   PACKED_MEMBER( uint32_t,     ulTimeSinceLastRF             );    // Sequence time stamp since last RF pulse           4
%   PACKED_MEMBER( uint16_t,     ushKSpaceCentreLineNo         );    // number of K-space centre line                     2
%   PACKED_MEMBER( uint16_t,     ushKSpaceCentrePartitionNo    );    // number of K-space centre partition                2
%   PACKED_MEMBER( uint16_t,     aushIceProgramPara[MDH_NUMBEROFICEPROGRAMPARA] ); // free parameter for IceProgram   8   =88
%   PACKED_MEMBER( uint16_t,     aushFreePara[MDH_FREEHDRPARA] );    // free parameter                          4 * 2 =   8
%   sSliceData sSD;                                                                                                  // Slice Data                                       28   =124
%   PACKED_MEMBER( uint16_t,       ushChannelId                  );    // channel Id must be the last parameter             2
%   PACKED_MEMBER( uint16_t,       ushPTABPosNeg                 );    // negative, absolute PTAB position in [0.1 mm]      2
%                                                                    // (automatically set by PCI_TX firmware)
% } sMDH;                                                            // total length: 32 * 32 Bit (128 Byte)            128
%
% sLC: ICE loopcounter: [Line Acq Slice Partition Echo Phase Rep Set Seg Ida Idb Idc Idd Ide]
%
% EvalInfoMask
% /*--------------------------------------------------------------------------*/
% /*  Definition of EvalInfoMask:                                             */
% /*--------------------------------------------------------------------------*/
% const MdhBitField MDH_ACQEND            ((unsigned long)0);
% const MdhBitField MDH_RTFEEDBACK        (1);
% const MdhBitField MDH_HPFEEDBACK        (2);
% const MdhBitField MDH_ONLINE            (3);
% const MdhBitField MDH_OFFLINE           (4);
% const MdhBitField MDH_SYNCDATA          (5);       // readout contains synchroneous data
% const MdhBitField MDH_LASTSCANINCONCAT  (8);       // Flag for last scan in concatination
% 
% const MdhBitField MDH_RAWDATACORRECTION (10);      // Correct the rawadata with the rawdata correction factor
% const MdhBitField MDH_LASTSCANINMEAS    (11);      // Flag for last scan in measurement
% const MdhBitField MDH_SCANSCALEFACTOR   (12);      // Flag for scan specific additional scale factor
% const MdhBitField MDH_2NDHADAMARPULSE   (13);      // 2nd RF exitation of HADAMAR
% const MdhBitField MDH_REFPHASESTABSCAN  (14);      // reference phase stabilization scan         
% const MdhBitField MDH_PHASESTABSCAN     (15);      // phase stabilization scan
% const MdhBitField MDH_D3FFT             (16);      // execute 3D FFT         
% const MdhBitField MDH_SIGNREV           (17);      // sign reversal
% const MdhBitField MDH_PHASEFFT          (18);      // execute phase fft     
% const MdhBitField MDH_SWAPPED           (19);      // swapped phase/readout direction
% const MdhBitField MDH_POSTSHAREDLINE    (20);      // shared line               
% const MdhBitField MDH_PHASCOR           (21);      // phase correction data    
% const MdhBitField MDH_PATREFSCAN        (22);      // additonal scan for PAT reference line/partition
% const MdhBitField MDH_PATREFANDIMASCAN  (23);      // additonal scan for PAT reference line/partition that is also used as image scan
% const MdhBitField MDH_REFLECT           (24);      // reflect line              
% const MdhBitField MDH_NOISEADJSCAN      (25);      // noise adjust scan --> Not used in NUM4        
% const MdhBitField MDH_SHARENOW          (26);      // all lines are acquired from the actual and previous e.g. phases
% const MdhBitField MDH_LASTMEASUREDLINE  (27);      // indicates that the current line is the last measured line of all succeeding e.g. phases
% const MdhBitField MDH_FIRSTSCANINSLICE  (28);      // indicates first scan in slice (needed for time stamps)
% const MdhBitField MDH_LASTSCANINSLICE   (29);      // indicates  last scan in slice (needed for time stamps)
% const MdhBitField MDH_TREFFECTIVEBEGIN  (30);      // indicates the begin time stamp for TReff (triggered measurement)
% const MdhBitField MDH_TREFFECTIVEEND    (31);      // indicates the   end time stamp for TReff (triggered measurement)
% const MdhBitField MDH_MDS_REF_POSITION  (32);      // indicates the reference position for move during scan images (must be set once per slice/partition in MDS mode)
% const MdhBitField MDH_SLC_AVERAGED      (33);      // indicates avveraged slice for slice partial averaging scheme
% 
% const MdhBitField MDH_FIRST_SCAN_IN_BLADE       (40);  // Marks the first line of a blade
% const MdhBitField MDH_LAST_SCAN_IN_BLADE        (41);  // Marks the last line of a blade
% const MdhBitField MDH_LAST_BLADE_IN_TR          (42);  // Set for all lines of the last BLADE in each TR interval
% 
% const MdhBitField MDH_RETRO_LASTPHASE           (45);  // Marks the last phase in a heartbeat
% const MdhBitField MDH_RETRO_ENDOFMEAS           (46);  // Marks an ADC at the end of the measurement
% const MdhBitField MDH_RETRO_REPEATTHISHEARTBEAT (47);  // Repeat the current heartbeat when this bit is found
% const MdhBitField MDH_RETRO_REPEATPREVHEARTBEAT (48);  // Repeat the previous heartbeat when this bit is found
% const MdhBitField MDH_RETRO_ABORTSCANNOW        (49);  // Just abort everything
% const MdhBitField MDH_RETRO_LASTHEARTBEAT       (50);  // This adc is from the last heartbeat (a dummy)
% const MdhBitField MDH_RETRO_DUMMYSCAN           (51);  // This adc is just a dummy scan, throw it away
% const MdhBitField MDH_RETRO_ARRDETDISABLED      (52);  // Disable all arrhythmia detection when this bit is found
% 
% //-----------------------------------------------------------------------------
% // Definition of EvalInfoMask for COP:
% //-----------------------------------------------------------------------------
% #define   MDH_COP_ACQEND            (0x00000001)
% #define MDH_COP_RTFEEDBACK        (0x00000002)
% #define MDH_COP_HPFEEDBACK        (0x00000004)
% #define   MDH_COP_ONLINE            (0x00000008)
% #define   MDH_COP_OFFLINE           (0x00000010)
 
% tic

measurementNum = size(measData, 1);
kspace = cell(measurementNum, 1);
Noise = cell(measurementNum, 1);
ref = cell(measurementNum, 1);

% [Data, asc, prot] =Read_RawData(fname);
for mi = 1 : measurementNum
        
    measData{mi}.prot = char(measData{mi}.prot); % header in text mode
    
    % get the FFT scale factor
    scaleFactor=findFFTScaleFactor_VD11( measData{mi}.prot);
    
    if scaleFactor==-1
        scaleFactor=ones(measData{mi}.asc(1).ushUsedChannels,1);
    end
    % if we did not find scale factor, 1 is used for all channels
   
    S0 = size(measData{mi}.asc);
    numOfLines = S0(2);
       
    % The actual channel id can be arbitrary between 0 and max.
    % available coil number. For recon, we can squeeze them to start
    % from 1, and consecutively to the actual used channel number
    measData{mi}.asc = reorderChannelId(measData{mi}.asc);

    % find the noise lines and seperate ref lines if any
    numOfNoiseLines = 0; 
    NoiseLines = zeros(numOfLines,1, 'int32');
    sLCNoise = zeros(numOfLines, 11, 'int32');

    numOfRefLines = 0; 
    RefLines = zeros(numOfLines,1, 'int32');
    sLCRef = zeros(numOfLines, 11, 'int32');

    SeparateRefLines = zeros(numOfLines,1, 'uint8');
    
    numOfDataLines = 0; 
    DataLines = zeros(numOfLines,1, 'int32');
    sLC = zeros(numOfLines, 11, 'int32');

    kspaceCentreLineNo = zeros(numOfLines,1, 'int32');

    for i=1:numOfLines
        
        if ( IsAcqEnd(measData{mi}.asc(i).aulEvalInfoMask(1)))
            continue;
        end
        
        kspaceCentreLineNo(i) = measData{mi}.asc(i).ushKSpaceCentreLineNo;

        [bNoiseLine, bSeparateRef, bRef, bOnline] = parseEvalInfoMask(measData{mi}.asc(i).aulEvalInfoMask(1));

        if ( bNoiseLine  )
            numOfNoiseLines = numOfNoiseLines + 1;
            NoiseLines(numOfNoiseLines) = i;
            sLCNoise(numOfNoiseLines,:) = fillSLC(measData{mi}.asc(i));
            continue;
        end

        if ( bRef )
            numOfRefLines = numOfRefLines + 1;
            RefLines(numOfRefLines) = i;
            sLCRef(numOfRefLines,:) = fillSLC(measData{mi}.asc(i));
            if (~flagIncludeSeparateRef && bSeparateRef) % bRef is true as long as bSeparateRef is true
                SeparateRefLines(i) = 1;
                continue;
            end
        end 

        if (bOnline)
            numOfDataLines = numOfDataLines + 1;
            DataLines(numOfDataLines) = i;
            sLC(numOfDataLines,:) = fillSLC(measData{mi}.asc(i));
        end
    end
    
    %% temporaly remove ref, as it is included in the data
    % change to "if 0" to disable this
    if 1
        numOfRefLines=0;
    end
    
    
    sLCNoise = sLCNoise(1:numOfNoiseLines, :); NoiseLines = NoiseLines(1:numOfNoiseLines);
    sLCRef = sLCRef(1:numOfRefLines, :); RefLines = RefLines(1:numOfRefLines);
    sLCRef(:,11) = 1; % Put different segments into one target k-space
    sLC = sLC(1:numOfDataLines, :); DataLines = DataLines(1:numOfDataLines);
    maxKSpaceLineNo = 2*max(kspaceCentreLineNo);
    sLC(:,11) = 1;         % Put different segments into one target k-space

    % find intended image width and height below
    [PElineOut, FElineOut] = findImageSizeFromHeader(measData{mi}.prot);

    % kspace size in  dim = [maxCol, maxLine, maxCha, maxAcq, maxSlice, maxPar,
    % maxEcho, maxPhs, maxRep, maxSet, maxSeg];
    kspaceDim = getOrgKSpaceSize(sLC, maxKSpaceLineNo);
    maxFElines = kspaceDim(1);
    FEbuffer = zeros(maxFElines,1);
    if flagUndoOS
        kspaceDim(1) = FElineOut;
        ROIstartIndex = (maxFElines - FElineOut)/2+1;
    end
    % Always use the intended PEline number to form the final k-space?
    % The difference between the computed kspaceDim(2) and intended PEline
    % number shouldn't be too large (give a warning if this happens)
%     if abs(kspaceDim(2) - PElineOut) > PElineOut/4
%         disp('Warning: the intended output imageLines is very different from the computed maxKSpaceLineNo');
%     end
    PElineOut= max(PElineOut,kspaceDim(2));
    kspaceDim(2) =PElineOut;
    outDCindex = ceil((PElineOut-1)/2); % Will be used to prevent phase shift. No impact on magnitude images
                                          % remove + 1 here by Jun Liu,
                                          % October 1, 2011
    kspace{mi}=complex(zeros(kspaceDim, 'single'), zeros(kspaceDim, 'single') );

    % [Col Line Cha Acq Slice Partition Echo Phase Rep Set Seg]
    % 11 dimension array

    % noise data if any
    if ( numOfNoiseLines > 0 )
        maxColNoise = max(sLCNoise(:, 1));
        Noise{mi} = zeros(maxColNoise, numOfNoiseLines);
    else
        Noise{mi} = [];
    end

    % seperate reference lines
    if ( numOfRefLines > 0 )
        ref{mi} =complex(zeros(kspaceDim, 'single'), zeros(kspaceDim, 'single') );
    else
        ref{mi} = [];
    end

    % put the data into the array

    noiseLineCount = 1;
    totalLineCount = 1;
    for i=1:numOfLines
        
        if ( mod(i, 1e2) == 0 )
            disp([num2str(i) ' lines parsed ... ']);
        end
        %fprintf('\n %d out of %d lines parsed...',i,numOfLines);
        %toc;

        if ( IsAcqEnd(measData{mi}.asc(i).aulEvalInfoMask(1)) )
            continue;
        end

        % the totalLineCount should be dependent on i only
        totalLineCount=(i-1) *kspaceDim(3) +1; % kspaceDim(3) is the number of coils
                
        if ( ~isempty(find(i==NoiseLines)) )
            FEbuffer(:) = 0;
            N = measData{mi}.asc(i).ushSamplesInScan;
            aZ = addPrePostZeros(measData{mi}.asc(i));
            % Take care of different channels here
            for ci = 1 : measData{mi}.asc(i).ushUsedChannels
                if ( aZ == 0 || aZ == 2)
                    FEbuffer(1:N) = measData{mi}.Data(1:N, totalLineCount);
                else %if ( aZ == 1 ) % pre zeros
                    FEbuffer(end-N+1:end) = measData{mi}.Data(1:N, totalLineCount);
                end
                if flagUndoOS % undo oversampling here
                    temp = fftshift(ifft(fftshift(FEbuffer)));
                    temp = temp(ROIstartIndex : ROIstartIndex+FElineOut-1);
                    temp = fftshift(fft(fftshift(temp)));                
                    Noise{mi}(1:FElineOut, noiseLineCount) = temp* scaleFactor(ci);
                else %if ( aZ == 1 ) % pre zeros
                    Noise{mi}(1:N, noiseLineCount) = FEbuffer* scaleFactor(ci);
                end
                noiseLineCount = noiseLineCount + 1;
                totalLineCount = totalLineCount + 1;
            end
            continue;
        end

        % the totalLineCount should be dependent on i only
        totalLineCount=(i-1) *kspaceDim(3) +1; % kspaceDim(3) is the number of coils
        
        r = find(i==RefLines);
        if ( ~isempty(r)  )
            FEbuffer(:) = 0;
            N = measData{mi}.asc(i).ushSamplesInScan;
            aZ = addPrePostZeros(measData{mi}.asc(i));
            PEindex = sLCRef(r,2);
            DCshift = outDCindex - kspaceCentreLineNo(i);
            outPEindex = PEindex + DCshift; % for correcting phase in the output k-space
            if ( outPEindex <= PElineOut && outPEindex > 0 )
            % Take care of different channels here
                for ci = 1 : measData{mi}.asc(i).ushUsedChannels
                    if ( aZ == 0 || aZ == 2 )
                        FEbuffer(1:N) = measData{mi}.Data(1:N, totalLineCount);
                    else % aZ == 1, pre-zeros
                        FEbuffer(end-N+1:end) = measData{mi}.Data(1:N, totalLineCount);
                    end
                    % adjust DC to the center of the output k-space
                    if flagUndoOS % undo oversampling here
                        temp = fftshift(ifft(fftshift(FEbuffer)));
                        temp = temp(ROIstartIndex : ROIstartIndex+FElineOut-1);
                        temp = fftshift(fft(fftshift(temp)));                        
                        ref{mi}(1:FElineOut, outPEindex, measData{mi}.asc(i).ushChannelId(ci), sLCRef(r,4), sLCRef(r,5), ...
                            sLCRef(r,6), sLCRef(r,7), sLCRef(r,8), sLCRef(r,9), sLCRef(r,10), sLCRef(r,11)) = temp * scaleFactor(ci);
                    else
                        ref{mi}(1:N, outPEindex, measData{mi}.asc(i).ushChannelId(ci), sLCRef(r,4), sLCRef(r,5), ...
                            sLCRef(r,6), sLCRef(r,7), sLCRef(r,8), sLCRef(r,9), sLCRef(r,10), sLCRef(r,11)) = FEbuffer* scaleFactor(ci);
                    end
                    totalLineCount = totalLineCount + 1;
                end
                if (~flagIncludeSeparateRef && SeparateRefLines(i))
                    continue;
                end
                totalLineCount = totalLineCount -  measData{mi}.asc(i).ushUsedChannels; % So that DataLines will use the same ref as DataLines
            end
        end

        % the totalLineCount should be dependent on i only
        totalLineCount=(i-1) *kspaceDim(3) +1; % kspaceDim(3) is the number of coils
        
        l = find(i==DataLines);
        if ( ~isempty(l) )
            FEbuffer(:) = 0;
            N = measData{mi}.asc(i).ushSamplesInScan;
            aZ = addPrePostZeros(measData{mi}.asc(i));
            PEindex = sLC(l,2);
            DCshift = outDCindex - kspaceCentreLineNo(i);
            outPEindex = PEindex + DCshift; % for correcting phase in the output k-space
            if ( outPEindex <= PElineOut && outPEindex > 0 )
                % Take care of different channels here
                for ci = 1 : measData{mi}.asc(i).ushUsedChannels
                    if ( aZ == 0 || aZ == 2) % no filling or post zeros
                        FEbuffer(1:N) = measData{mi}.Data(1:N, totalLineCount);
                    else % aZ == 1, pre-zeros
                        FEbuffer(end-N+1:end) = measData{mi}.Data(1:N, totalLineCount);
                    end
                    % adjust DC to the center of the output k-space
                    if flagUndoOS % undo oversampling here
                        temp = fftshift(ifft(fftshift(FEbuffer)));
                        temp = temp(ROIstartIndex : ROIstartIndex+FElineOut-1);
                        temp = fftshift(fft(fftshift(temp)));
                        kspace{mi}(1:FElineOut, outPEindex, measData{mi}.asc(i).ushChannelId(ci), sLC(l,4), sLC(l,5), ...
                                sLC(l,6), sLC(l,7), sLC(l,8), sLC(l,9), sLC(l,10), sLC(l,11)) = temp* scaleFactor(ci);
                    else
                        kspace{mi}(1:N, outPEindex, measData{mi}.asc(i).ushChannelId(ci), sLC(l,4), sLC(l,5), ...
                                sLC(l,6), sLC(l,7), sLC(l,8), sLC(l,9), sLC(l,10), sLC(l,11)) = FEbuffer* scaleFactor(ci);
                    end
                    totalLineCount = totalLineCount + 1;                    
                end
            end % of ( sLC(l,2) <= maxKSpaceLineNo )
        end
    end
end

% toc

%% Used function below
    % --------------------------------------------
    % parse the Eval Info Mask
    function bValue = IsNoiseScan(evalInfoMask)
        % const MdhBitField MDH_NOISEADJSCAN      (25);      // noise adjust scan --> Not used in NUM4
        p = dec2bin(evalInfoMask(1));
        len = numel(p);
        if ( len >= 26 )
            if ( p(end-25) == '1' )
                bValue = 1;
            else
                bValue = 0;
            end
        else
            bValue = 0;
        end
    end
 
    % --------------------------------------------
    % the seperate reference lines
    function bValue = IsSeperateRef(evalInfoMask)
        % const MdhBitField MDH_PATREFSCAN        (22);      // additonal scan for PAT reference line/partition
        p = dec2bin(evalInfoMask(1));
        len = numel(p);
        if ( len >= 23 )
            if ( p(end-22) == '1' )
                bValue = 1;
            else
                bValue = 0;
            end
        else
            bValue = 0;
        end
    end
    % --------------------------------------------
    % parse the Eval Info Mask
    function [bNoiseLine, bSeparateRef, bRef, bOnline] = parseEvalInfoMask(evalInfoMask)
        % const MdhBitField MDH_NOISEADJSCAN      (25);      // noise adjust scan --> Not used in NUM4
        % const MdhBitField MDH_PATREFSCAN        (22);      // additonal scan for PAT reference line/partition
        % const MdhBitField MDH_PATREFANDIMASCAN  (23);      // additonal scan for PAT reference line/partition that is also used as image scan
        infoMask = ['00000000000000000000000000000000'];
        p = dec2bin(evalInfoMask(1));
        len = numel(p);
        
        
        
        if ( len < 32 )
            infoMask(end-len+1:end) = p;
        else
            infoMask = p;
        end
 
        bNoiseLine = 0;
        if ( infoMask(end-25) == '1' )
            bNoiseLine = 1;
        end
        bSeparateRef = 0;
        if ( infoMask(end-22)=='1' && infoMask(end-23)~='1' )
            bSeparateRef = 1;
        end
        
        bRef = 0;
        if ( infoMask(end-22)=='1' || infoMask(end-23)=='1' )
            bRef = 1;
        end
        bOnline=0;
        if (infoMask(end-3) == '1')
            bOnline = 1;
        end
    end
    % --------------------------------------------
    % the acq end line
    function bValue = IsAcqEnd(evalInfoMask)
        % const MdhBitField MDH_ACQEND            ((unsigned long)0);
        bValue = (evalInfoMask(1)==1);
    end
    % --------------------------------------------

    % get kSpace size
    function dim = getOrgKSpaceSize(sLC, maxKSpaceLineNo)
        maxCol = max(sLC(:,1)); % starting from 1
        % maxLine = min(max(sLC(:,2)), maxKSpaceLineNo);
        maxLine = max(max(sLC(:,2)), maxKSpaceLineNo);
        maxCha = max(sLC(:,3));
        maxAcq = max(sLC(:,4));
        maxSlice = max(sLC(:,5));
        maxPar = max(sLC(:,6));
        maxEcho = max(sLC(:,7));
        maxPhs = max(sLC(:,8));
        maxRep = max(sLC(:,9));
        maxSet = max(sLC(:,10));
        maxSeg = max(sLC(:,11));
 
        % combine the data dimension size
        dim = [maxCol, maxLine, maxCha, maxAcq, maxSlice, maxPar, maxEcho, maxPhs, maxRep, maxSet, maxSeg];
    end

    % --------------------------------------------
    % fill the sLC
    function aSLC = fillSLC(aAsc)
        
        if ( aAsc.ushKSpaceCentreColumn == 0 )
            aSLC(1) = aAsc.ushSamplesInScan;
        else
            if ( 2*aAsc.ushKSpaceCentreColumn >= aAsc.ushSamplesInScan )
                aSLC(1) = 2*aAsc.ushKSpaceCentreColumn;
            else
                aSLC(1) = 2*(aAsc.ushSamplesInScan - aAsc.ushKSpaceCentreColumn);
            end
        end
%         aSLC(1) = max([2*aAsc.ushKSpaceCentreColumn aAsc.ushSamplesInScan]);
        aSLC(2) = aAsc.sLC(1)+1; % Line
%         aSLC(3) = aAsc.ushChannelId+1; % channel
        aSLC(3) = aAsc.ushUsedChannels;
        aSLC(4:11) = aAsc.sLC(2:9)+1; % other dimensions
    end
    % --------------------------------------------
    % add pre or post zeros
    function aZ = addPrePostZeros(aAsc)
        % aZ = 1 : pre zeros
        % aZ = 2 : post zeros
        % aZ = 0 : no zeros
        if ( 2*aAsc.ushKSpaceCentreColumn == aAsc.ushSamplesInScan )
            aZ = 0;
            return;
        end
        
        if ( 2*aAsc.ushKSpaceCentreColumn < aAsc.ushSamplesInScan )
            aZ = 1;
            return;
        end
        
        if ( 2*aAsc.ushKSpaceCentreColumn > aAsc.ushSamplesInScan )
            aZ = 2;
            return;
        end
        
    end



        function flagIsValidLine = checkIfValidLine(asci)
            if asci.ushSamplesInScan>0
                flagIsValidLine = true;
            else
                flagIsValidLine = false;
            end
        end
 
end
 
