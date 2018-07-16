function [out] = read_vd_datfile(datfilename, measurement, flag);
% function [out] = read_vd_datfile(datfilename, measurement, flag);
%
% function to read siemens raw data file(VD version)
%
% input:
%     datfilename  name of raw data file
%     measurement  1,2,3  index that specified which multiraid dataset
%                  (0 indicates to read the last measurement)
%     flag: 0, 1(def)  0: reads only headers, PMU data, readoutindicators & 1 reads raw data as well (Slower)
% output:
%
%     data(encode, freq, acq, slice, partition, echo, phase,repetition, dataset, segment, Ida, Idb, Idc, Idd, Ide, coil)
%     readoutindicator is a matrix with a row vector for each readout containing:
%     readoutindicator(i,:)= [ ...
%         Line,...       % line index                   */
%         Acquisition,...% acquisition index            */
%         Slice,...      % slice index                  */
%         Partition,...  % partition index              */
%         Echo,...       % echo index                   */
%         Phase,...      % phase index                  */
%         Repetition,... % measurement repeat index     */
%         Set,...        % set index                    */
%         Seg,...        % segment index  (for TSE)     */
%         EvalInfoMask,...%
%         TimeStamp,...       % time stamp [2.5 ms ticks since 00:00]
%         PMUTimeStamp,...    % PMU time stamp [2.5 ms ticks since last trigger]
%         SamplesInScan,...  % # of samples acquired in scan
%         Ida,...        % IceDimension a index         */
%         Idb,...        % IceDimension b index         */
%         Idc,...        % IceDimension c index         */
%         Idd,...        % IceDimension d index         */
%         Ide];          % IceDimension e index         */
%
%              for example: plot(readoutindicator(:,1),'.') plots the phase encode order
%                           plot(2.5e-3*readoutindicator(:,11),'.') plots the time stamp of each readout
    
%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  NIH NHLBI                          *
%     ***************************************

% set defaults
if nargin < 2
    measurement = 0; % will choose last measurement
end
if nargin < 3
    flag = 1; % read data by default
end

verbose = 0;
[tmp,name] = fileparts(datfilename);
out.datfilename = name;
skip_sync_data = 0;
[MrParcRaidFileHeader, MrParcRaidFileEntry, headers, protocol]=read_vd_datfile_headers([datfilename]);
out.headers = headers;
out.protocol = protocol;

Nmeasurements = length(MrParcRaidFileEntry);
if measurement == 0
    measurement = Nmeasurements;
end


fid=fopen(datfilename,'r','n');

fseek(fid,MrParcRaidFileEntry(measurement).off,-1);% skip to file offset for beginning of measurement headers
dma_len = fread(fid,1,'uint32'); % read the length of the measurement header including alignment.

fseek(fid, dma_len - 4 ,0); % jump to beginning of "Measurement"

MDH_ACQEND = 0;

i = 1; j = 1; k = 1;
while MDH_ACQEND == 0
    % read sScanHeader (192 bytes)
    if verbose
        disp(['k: ',num2str(k)])
    end
    CurrentScanHeader = read_n4_header(fid);

    if verbose
        disp(['MeasUID: ',num2str(CurrentScanHeader.lMeasUID)])
        disp(['CurrentScanHeader.aulEvalInfoMask(1): ',num2str(CurrentScanHeader.aulEvalInfoMask(1))])
        disp(['CurrentScanHeader.aulEvalInfoMask(2): ',num2str(CurrentScanHeader.aulEvalInfoMask(2))])
    end
    
%     disp(num2str(ScanHeader(i).aulEvalInfoMask))
    if numel(CurrentScanHeader.aulEvalInfoMask)==0; ScanHeader = ScanHeader(1:end-1); break; end
    MDH_ACQEND  = bitand(CurrentScanHeader.aulEvalInfoMask(1), 2^0);
    MDH_SYNCDATA  = bitand(CurrentScanHeader.aulEvalInfoMask(1), 2^5);
    MDH_LASTSCANINSLICE = bitand(CurrentScanHeader.aulEvalInfoMask(1), 2^29);
    tmp = [];
   
    if verbose
        disp(['MDH_ACQEND: ',num2str(MDH_ACQEND)])
    end
    
    if MDH_ACQEND == 0
    if MDH_SYNCDATA == 0 
        if flag == -1;
            out.ScanHeader = CurrentScanHeader;
            return
        end
        ScanHeader(i) = CurrentScanHeader;
        ushUsedChannels = CurrentScanHeader.ushUsedChannels;
        for c = 1:ushUsedChannels % read ChannelData
            % read sChannelHeader (32 bytes)
            ulTypeAndChannelLength = fread(fid,1,'uint32');
            lMeasUID = fread(fid,1,'uint32');
            ulScanCounter = fread(fid,1,'uint32');
            ulReserved1 = fread(fid,1,'uint32');
            ulSequenceTime = fread(fid,1,'uint32');
            ulUnused2 = fread(fid,1,'uint32');
            ulChannelID = fread(fid,1,'uint16');
            ulUnused3 = fread(fid,1,'uint16');
            ulCRC = fread(fid,1,'uint32');
            % read raw data
            tmp(:,c) = fread(fid,2*CurrentScanHeader.ushSamplesInScan,'single');
        end
        if ushUsedChannels > 0
            rawdata{i} = complex(tmp(1:2:end,:),tmp(2:2:end,:));
            if verbose; disp([num2str(i),'  ',num2str(CurrentScanHeader.ushSamplesInScan)]); end
            i = i + 1;
        end
    else % is SYNC data (containing PMU signals)
            PMU.ScanHeader(j) = CurrentScanHeader;
            ulFlagsAndDMALength = CurrentScanHeader.ulFlagsAndDMALength;
            % read SeqDataHeader (60 bytes)
            pos1 = ftell(fid);
            PacketSize = fread(fid,1,'uint32'); % size includes header
            PMU.PackedID{j} = char(fread(fid,52,'char')');
            if verbose
                disp(['PackedID:   ',PMU.PackedID{j}])
                disp(['PacketSize: ',num2str(PacketSize)])
            end
            
            if isempty(findstr(PMU.PackedID{j},'PMU')) % then not PMU data or PMULearnedphase
                skip_sync_data = 1;
                clear PMU
            end
            
            SwappedFlag = fread(fid,1,'uint32');
            if skip_sync_data == 1
                dummy_bytes = ceil((PacketSize+60)/32)*32 - PacketSize - 60; 
                fseek(fid, PacketSize + dummy_bytes  ,0); % jump to beginning of next scan
            else

                % read PMU data
                % end indicated when MagicEnd  = PMU_MAGIC_END = 0x01FF0000 (33488896)
                % PMU_MAGIC_ECG1 = 0x01010000 (16842752)
                % PMU_MAGIC_ECG2 = 0x01020000 (16908288)
                % PMU_MAGIC_ECG3 = 0x01030000 (16973824)
                % PMU_MAGIC_ECG4 = 0x01040000 (17039360)
                % PMU_MAGIC_PULS = 0x01050000 (17104896)
                % PMU_MAGIC_RESP = 0x01060000 (17170432)
                % PMU_MAGIC_EXT1 = 0x01070000
                % PMU_MAGIC_EXT2 = 0x01080000
                
                PMU.TimeStamp(j,:) = fread(fid,2,'uint32');                
                PMU.PackerNr(j) = fread(fid,1,'uint32');
                if verbose; disp(['PackerNr: ',num2str(PMU.PackerNr(j))]);end
                Duration = fread(fid,1,'uint32');
                Magic = 0;
                while Magic ~= 33488896;
                    Magic = fread(fid,1,'uint32');
                    if Magic ~= 33488896; Period = fread(fid,1,'uint32'); end
                    switch Magic
                        case 16842752 % ECG1
                            PMU.Data{1,j} = fread(fid,Duration/Period,'uint32');
                            PMU.ECGDataPeriod{j} = Period;
                            PMU.ECGDataDuration{j} = Duration;
                        case 16908288 % ECG2
                            PMU.Data{2,j} = fread(fid,Duration/Period,'uint32');
                        case 16973824 % ECG3
                            PMU.Data{3,j} = fread(fid,Duration/Period,'uint32');
                        case 17039360 % ECG4
                            PMU.Data{4,j} = fread(fid,Duration/Period,'uint32');
                        case 17104896 % PULS
                            PMU.Data{5,j} = fread(fid,Duration/Period,'uint32');
                        case 17170432 % RESP
                            PMU.Data{6,j} = fread(fid,Duration/Period,'uint32');                            
                            PMU.RESPDataPeriod{j} = Period;
                            PMU.RESPDataDuration{j} = Duration;
                    end
            
                end
                pos2 = ftell(fid);
                j=j+1;
                dummy_bytes = ceil((pos2-pos1)/32)*32 - (pos2-pos1); 
                fseek(fid, dummy_bytes  ,0); % jump to beginning of next scan
                if verbose
                    disp(['ulFlagsAndDMALength: ',num2str(ulFlagsAndDMALength)])
                    disp(['dummy_bytes: ',num2str(dummy_bytes)])
                    disp(['calculated: ',num2str(192+PacketSize+60+dummy_bytes)])
                end
            end
        end
    k = k + 1;
    if verbose; disp(' '); end
    else
        if verbose; disp(['MDH_ACQEND=1']); end
        if verbose; disp(' '); end
    end
%     if MDH_LASTSCANINSLICE; break; end
end
     
fclose(fid);


% d = complex(tmp(1:2:end,:),tmp(2:2:end,:));
for i = 1:length(ScanHeader)
    readoutindicator(i,:)= [ ...
        ScanHeader(i).sLC.ushLine,...       % line index                   */
        ScanHeader(i).sLC.ushAcquisition,...% acquisition index            */
        ScanHeader(i).sLC.ushSlice,...      % slice index                  */
        ScanHeader(i).sLC.ushPartition,...  % partition index              */
        ScanHeader(i).sLC.ushEcho,...       % echo index                   */
        ScanHeader(i).sLC.ushPhase,...      % phase index                  */
        ScanHeader(i).sLC.ushRepetition,... % measurement repeat index     */
        ScanHeader(i).sLC.ushSet,...        % set index                    */
        ScanHeader(i).sLC.ushSeg,...        % segment index  (for TSE)     */
        ScanHeader(i).aulEvalInfoMask(1),...%
        ScanHeader(i).ulTimeStamp,...       % time stamp [2.5 ms ticks since 00:00]
        ScanHeader(i).ulPMUTimeStamp,...    % PMU time stamp [2.5 ms ticks since last trigger]
        ScanHeader(i).ushSamplesInScan,...  % # of samples acquired in scan
        ScanHeader(i).sLC.ushIda,...        % IceDimension a index         */
        ScanHeader(i).sLC.ushIdb,...        % IceDimension b index         */
        ScanHeader(i).sLC.ushIdc,...        % IceDimension c index         */
        ScanHeader(i).sLC.ushIdd,...        % IceDimension d index         */
        ScanHeader(i).sLC.ushIde,...        % IceDimension e index         */
        ScanHeader(i).ulTimeSinceLastRF];
end
out.ScanHeader = ScanHeader;
out.readoutindicator = readoutindicator;

aulEvalInfoMask = [ScanHeader(:).aulEvalInfoMask];
MDH_NOISEADJSCAN  = bitand(aulEvalInfoMask(1,:), 2^25);
MDH_PATREFSCAN  = bitand(aulEvalInfoMask(1,:), 2^22);
MDH_PATREFANDIMASCAN  = bitand(aulEvalInfoMask(1,:), 2^23);
MDH_REFLECT  = bitand(aulEvalInfoMask(1,:), 2^24);
MDH_PHASCOR  = bitand(aulEvalInfoMask(1,:), 2^21);
MDH_RTFEEDBACK  = bitand(aulEvalInfoMask(1,:), 2^1);
MDH_HPFEEDBACK  = bitand(aulEvalInfoMask(1,:), 2^2);
MDH_SYNCDATA  = bitand(aulEvalInfoMask(1,:), 2^5);
MDH_ACQEND  = bitand(aulEvalInfoMask(1,:), 2^0);
MDH_PPAEXTRAREFSCAN  = bitand(aulEvalInfoMask(1,:), 2^13);


noisescans = find(MDH_NOISEADJSCAN);
noise =[];
% for i = 1:length(noisescans)
%     noise(:,:,i) = rawdata{noisescans(i)}.';
% end
% noise = noise(:,:); % channels x samples
for i = 1:length(noisescans)
    noise = [noise,rawdata{noisescans(i)}.'];
end



switch protocol{end}.sPat.ucRefScanMode
    case 4 % separate reference lines
        patrefdata = [];
        patrefscans = find(MDH_PATREFSCAN);
        for i = 1:length(patrefscans)
            patrefdata(:,:,i) = rawdata{patrefscans(i)}.';
        end
        out.acs = permute(patrefdata,[3 2 1]);
end

ignore = MDH_NOISEADJSCAN | MDH_PHASCOR | MDH_RTFEEDBACK | MDH_HPFEEDBACK | MDH_SYNCDATA | MDH_ACQEND | MDH_PPAEXTRAREFSCAN;

ushLine            =  readoutindicator(:,1);
ushAcquisition     =  readoutindicator(:,2);
ushSlice           =  readoutindicator(:,3);
ushPartition       =  readoutindicator(:,4);
ushEcho            =  readoutindicator(:,5);
ushPhase           =  readoutindicator(:,6);
ushRepetition      =  readoutindicator(:,7);
ushSet             =  readoutindicator(:,8);
ushSeg             =  readoutindicator(:,9);
ushIda             =  readoutindicator(:,14)+1;
ushIdb             =  readoutindicator(:,15)+1;
ushIdc             =  readoutindicator(:,16)+1;
ushIdd             =  readoutindicator(:,17)+1;
ushIde             =  readoutindicator(:,18)+1;


ushSamplesInScan = readoutindicator(:,13)';
ushSamplesInScan_unique = unique(ushSamplesInScan(~ignore));

% nsamples = ushSamplesInScan_unique(1);
% imagescans(1) = find(~ignore & ushSamplesInScan == nsamples);
%skip = length(protocol{1}.sCoilSelectMeas.aRxCoilSelectData_); % used for interleaving of different sets of coils
skip = 1;

if flag == 1
for num = 1:length(ushSamplesInScan_unique)
    nsamples = ushSamplesInScan_unique(num);
    
    imagescans = find(~ignore & ushSamplesInScan == nsamples);

    npe            =  max(ushLine(imagescans));
    nacqs          =  max(ushAcquisition(imagescans));
    nslices        =  max(ushSlice(imagescans));
    npartitions    =  max(ushPartition(imagescans));
    necho          =  max(ushEcho(imagescans));
    nphases        =  max(ushPhase(imagescans));
    nrepetitions   =  max(ushRepetition(imagescans));
    ndatasets      =  max(ushSet(imagescans));
    nsegs          =  max(ushSeg(imagescans));
    nIda           =  max(ushIda(imagescans));
    nIdb           =  max(ushIdb(imagescans));
    nIdc           =  max(ushIdc(imagescans));
    nIdd           =  max(ushIdd(imagescans));
    nIde           =  max(ushIde(imagescans));


    [ncols, ncoils] = size(rawdata{imagescans(1)});
    data = zeros(ncols,ncoils,npe,nacqs,nslices,npartitions,necho,nphases,...
        nrepetitions,ndatasets,nsegs,nIda,nIdb,nIdc,nIdd,nIde, 'single');
    for i = 1:skip:length(imagescans)
        encode = ushLine(imagescans(i));
        acq = ushAcquisition(imagescans(i));
        slice = ushSlice(imagescans(i));
        partition = ushPartition(imagescans(i));
        echo = ushEcho(imagescans(i));
        phase = ushPhase(imagescans(i));
        repetition = ushRepetition(imagescans(i));
        dataset = ushSet(imagescans(i));
        segment = ushSeg(imagescans(i));
        Ida = ushIda(imagescans(i));
        Idb = ushIdb(imagescans(i));
        Idc = ushIdc(imagescans(i));
        Idd = ushIdd(imagescans(i));
        Ide = ushIde(imagescans(i));

        data(:,:,encode,acq,slice,partition,echo,phase,repetition,dataset,segment,Ida,Idb,Idc,Idd,Ide) = ...
            rawdata{imagescans(i)};
    end
    out.data{num} = permute(data,[3 1 4:ndims(data), 2]);
end

out.noise = noise;
% out.patrefdata = patrefdata;
else % flag ~=1
    out.data=[];
    out.noise=[];
    for num = 1:length(ushSamplesInScan_unique)
        nsamples = ushSamplesInScan_unique(num);
        imagescans = find(~ignore & ushSamplesInScan == nsamples);
    end
end
    
if exist('PMU')
    out.PMUraw = PMU;

    % bit 15 : table is moving
    % bit 14 : ECG trigger
    % bit 13 : Pulse trigger
    % bit 12 : Resp. trigger
    % bit 11 : EXT1 trigger
    % bit 10 : EXT2 trigger
    %     static const uint32_t TABLE_MOVING_MARK = 0x01 << 31;
    %     static const uint32_t ECG_TRIGGER_MARK  = 0x01 << 30;
    %     static const uint32_t PULS_TRIGGER_MARK = 0x01 << 29;
    %     static const uint32_t RESP_TRIGGER_MARK = 0x01 << 28;
    %     static const uint32_t EXT1_TRIGGER_MARK = 0x01 << 27;
    %     static const uint32_t EXT2_TRIGGER_MARK = 0x01 << 26;
    tmp = []; tmp_time = [];
    for i1 = 1:size(PMU.Data,2);
        if ~isempty(findstr(PMU.PackedID{i1},'PMUData'));
            tmp = [tmp,PMU.Data{1,i1}'];
            len = length(PMU.Data{1,i1});
            pmu_timestamp = PMU.TimeStamp(i1,2); % lower 32 bit word time stamp for j-th packet
            tmp_time = [tmp_time, pmu_timestamp*2.5e-3+ 2.5e-3*[0:len-1]]; % ECG samples are in 0.1ms increments
            % from pmu timestamp which is in 2.5e-3 increments from 00:00:00
        end;
    end;
    ecg1 = bitand(tmp,32767);% lower 16 bits
    trigger1 = bitand(tmp, 2^30);
    mask1 = bitand(tmp, hex2dec('FFFF0000'))/2^16; % upper 16 bits
    ecg1_time = tmp_time;

    tmp = [];for i1 = 1:size(PMU.Data,2); if ~isempty(findstr(PMU.PackedID{i1},'PMUData'));tmp = [tmp,PMU.Data{2,i1}'];end;end;
    ecg2 = bitand(tmp,32767);
    trigger2 = bitand(tmp, 2^30);

    tmp = [];for i1 = 1:size(PMU.Data,2); if ~isempty(findstr(PMU.PackedID{i1},'PMUData'));tmp = [tmp,PMU.Data{3,i1}'];end;end;
    ecg3 = bitand(tmp,32767);
    trigger3 = bitand(tmp, 2^30);

    tmp = [];for i1 = 1:size(PMU.Data,2); if ~isempty(findstr(PMU.PackedID{i1},'PMUData'));tmp = [tmp,PMU.Data{4,i1}'];end;end;
    ecg4 = bitand(tmp,32767);
    trigger4 = bitand(tmp, 2^30);

    tmp = [];for i1 = 1:size(PMU.Data,2); if ~isempty(findstr(PMU.PackedID{i1},'PMUData'));tmp = [tmp,PMU.Data{5,i1}'];end;end;
    puls = bitand(tmp,32767);
    trigger_puls = bitand(tmp, 2^29);

    tmp = [];for i1 = 1:size(PMU.Data,2); if ~isempty(findstr(PMU.PackedID{i1},'PMUData'));tmp = [tmp,PMU.Data{6,i1}'];end;end;
    resp = bitand(tmp,32767);
    trigger_resp = bitand(tmp, 2^28);


    ecg1 = (ecg1-2048); % output in mV
    ecg2 = (ecg2-2048);
    ecg3 = (ecg3-2048);
    ecg4 = (ecg4-2048);

    out.PMU.ecg1=ecg1;
    out.PMU.trigger_ecg1 = trigger1;
    out.PMU.mask1 = mask1;
    out.PMU.ecg1_time = ecg1_time;
    out.PMU.ecg2=ecg2;
    out.PMU.trigger_ecg2 = trigger2;
    out.PMU.ecg3=ecg3;
    out.PMU.trigger_ecg3 = trigger3;
    out.PMU.ecg4=ecg4;
    out.PMU.trigger_ecg4 = trigger4;
    out.PMU.puls=puls;
    out.PMU.trigger_puls = trigger_puls;
    out.PMU.resp=resp;
    out.PMU.trigger_resp = trigger_resp;

    out.PMUtraining = PMU_training(PMU);


    start_mdh = imagescans(1);
    t = zeros(1,length(out.PMUraw.PackedID));
    for i = 1:length(out.PMUraw.PackedID);
        t(i) = ~isempty(findstr(out.PMUraw.PackedID{i},'PMUData'));
    end
    start_pmu = min(find(t==1));
    if verbose
        disp(['start mdh: ',num2str(start_mdh)])
        disp(['start pmu: ',num2str(start_pmu)])
    end
    start_mdh=1;
    dt1 = (readoutindicator(end,11)-readoutindicator(start_mdh,11))*2.5e-3; % 2.5 msec increments
    dt2 = (out.PMUraw.TimeStamp(end,end)-out.PMUraw.TimeStamp(start_pmu+1,end))*100e-6;% 100 microsec increments
    if verbose
        disp(['dt1 (mdh):', num2str(dt1)])
        disp(['dt2 (pmu):', num2str(dt2)])
    end
end


return


RTFEEDBACK_flag=mask.MDH_RTFEEDBACK~=0;
MDH_ACQEND_flag=mask.MDH_ACQEND~=0;
MDH_PATREFSCAN_FLAG=mask.MDH_PATREFSCAN~=0;
patrefscan_flag = ((bitand(ulEvalInfoMask, 2^22)~=0)& (patrefscanmode>=1) );
valid = ~or(RTFEEDBACK_flag, MDH_ACQEND_flag);

switch protocol.sPat.ucRefScanMode
    case {' 0x4', ' 0x04', ' 0x2', ' 0x02'}
        output.pat_reference_data=pat_reference_data;
end





return

function [header] = read_n4_header(fid);
% function [header] = read_n4_header(fid,version);
%
% function to read siemens raw data file header (version Numaris 4)
%
% input:
%     fid     file handle, e.g.   fid=fopen(rawfilename,'r','n');
% output:
%     header  structure containing MDH header fields.

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************
header.ulFlagsAndDMALength = fread(fid,1,'int32'); % ulFlagsAndDMALength
header.lMeasUID            = fread(fid,1,'int32'); % measurement user ID
header.ulScanCounter       = fread(fid,1,'uint32');  % scan counter [1...]
header.ulTimeStamp         = fread(fid,1,'uint32');  % time stamp [2.5 ms ticks since 00:00]
header.ulPMUTimeStamp      = fread(fid,1,'uint32'); % PMU time stamp [2.5 ms ticks since last trigger]
header.ushSystemType       = fread(fid,1,'uint16'); %
header.ushPTABPosDelay     = fread(fid,1,'uint16'); %
header.lPTABPosX           = fread(fid,1,'uint32'); %
header.lPTABPosY           = fread(fid,1,'uint32'); %
header.lPTABPosZ           = fread(fid,1,'uint32'); %
header.ulReserved1         = fread(fid,1,'uint32'); %
header.aulEvalInfoMask     = fread(fid,2,'uint32'); %
header.ushSamplesInScan    = fread(fid,1,'uint16');  % # of samples acquired in scan
header.ushUsedChannels     = fread(fid,1,'uint16'); % # of channels used in scan
% loop counters
header.sLC.ushLine         = 1+fread(fid,1,'uint16'); % line index                   */
header.sLC.ushAcquisition  = 1+fread(fid,1,'uint16'); % acquisition index            */
header.sLC.ushSlice        = 1+fread(fid,1,'uint16'); % slice index                  */
header.sLC.ushPartition    = 1+fread(fid,1,'uint16'); % partition index              */
header.sLC.ushEcho         = 1+fread(fid,1,'uint16'); % echo index                   */
header.sLC.ushPhase        = 1+fread(fid,1,'uint16'); % phase index                  */
header.sLC.ushRepetition   = 1+fread(fid,1,'uint16'); % measurement repeat index     */
header.sLC.ushSet          = 1+fread(fid,1,'uint16'); % set index                    */
header.sLC.ushSeg          = 1+fread(fid,1,'uint16'); % segment index  (for TSE)     */
% header.sLC.ushFree         = fread(fid,1,'ushort'); %  free loop counter           */   
header.sLC.ushIda =fread(fid,1,'ushort');  % IceDimension a index        */
header.sLC.ushIdb =fread(fid,1,'ushort');  % IceDimension b index        */
header.sLC.ushIdc =fread(fid,1,'ushort');  % IceDimension c index        */
header.sLC.ushIdd =fread(fid,1,'ushort');  % IceDimension d index        */
header.sLC.ushIde =fread(fid,1,'ushort');  % IceDimension e index        */
% cut-off values  
header.sCutOffData.ushPre  = fread(fid,1,'uint16'); % write ushPre zeros at line start */
header.sCutOffData.ushPost = fread(fid,1,'uint16'); % write ushPost zeros at line end  */
% header.ushKSpaceCentreColumn   = fread(fid,1,'uint32');  % centre of echo
header.ushKSpaceCentreColumn   = fread(fid,1,'uint16');  % centre of echo
% header.ushCoilSelect       = fread(fid,1,'uint32');  % Bit 0..3: CoilSelect
header.ushCoilSelect       = fread(fid,1,'uint16');  % Bit 0..3: CoilSelect
% header.fReadOutOffcentre       = fread(fid,1,'int32');  % ReadOut offcenter value
header.fReadOutOffcentre       = fread(fid,1,'single');  % ReadOut offcenter value
header.ulTimeSinceLastRF       = fread(fid,1,'int32');   % Sequence time stamp since last RF pulse           
header.ushKSpaceCentreLineNo   = fread(fid,1,'uint16');  % number of K-space centre line
header.ushKSpaceCentrePartitionNo = fread(fid,1,'uint16');  % number of K-space centre partition
% slice data (28 bytes)
header.sSD.sVector.flSag =  fread(fid,1,'single');%  slice position vector        */
header.sSD.sVector.flCor =  fread(fid,1,'single');%  slice position vector        */
header.sSD.sVector.flTra =  fread(fid,1,'single');%  slice position vector        */
header.sSD.aflQuaternion.a =  fread(fid,1,'single'); % rotation matrix as quaternion*/
header.sSD.aflQuaternion.b =  fread(fid,1,'single'); % rotation matrix as quaternion*/
header.sSD.aflQuaternion.c =  fread(fid,1,'single'); % rotation matrix as quaternion*/
header.sSD.aflQuaternion.d =  fread(fid,1,'single'); % rotation matrix as quaternion*/
% other params
header.aushIceProgramPara=fread(fid,24,'uint16');
% header.aushReservedPara=fread(fid,2,'uint16');
header.aushReservedPara=fread(fid,4,'uint16');
header.ushApplicationCounter1 = fread(fid,1,'uint16');
header.ushApplicationCounter2 = fread(fid,1,'uint16');
header.ulCRC                  = fread(fid,1,'uint32');

return

function readoutindicator_valid_data=check_valid_data(readoutindicator);
% function readoutindicator_valid_data=check_valid_data(readoutindicator);
%
% function to delete rows (readouts) from readoutindicator which correspond
% to noiseadjscan or reference scan readouts

% define ulEvalInfoMask bit masks (from mdh.h include file)
mdh_noiseadjscan=base2dec('02000000',16);% define MDH_NOISEADJSCAN  0x02000000 // noise adjust scan
mdh_extrarefscan=base2dec('00002000',16);% define MDH_PPAEXTRAREFSCAN  0x00002000 // additional scan for PPA reference line/partition
mdh_phasecorrection=base2dec('00200000',16); % define MDH_PHASCOR 0x00200000 // phase correction data    
mdh_acqend=base2dec('00000001',16); % define MDH_ACQEND (end of acquisition dummy readout)

noiseadjscan_flag=bitand(readoutindicator(:,10),mdh_noiseadjscan)>0;
extrarefscan_flag=bitand(readoutindicator(:,10),mdh_extrarefscan)>0;
phasecorrection_flag=bitand(readoutindicator(:,10),mdh_phasecorrection)>0;
acqend_flag=bitand(readoutindicator(:,10),mdh_acqend)>0;
 
valid_flag=(noiseadjscan_flag==0) & (extrarefscan_flag==0) & ...
    (phasecorrection_flag==0) & (acqend_flag==0);
% valid_flag=(noiseadjscan_flag==0) & ...
%     (phasecorrection_flag==0) & (acqend_flag==0);
validreadouts=find(valid_flag==1);
readoutindicator_valid_data=readoutindicator(validreadouts,:);

return

% const MdhBitField MDH_ACQEND            (0UL);       ///< last scan 
% const MdhBitField MDH_RTFEEDBACK        (1UL);       ///< Realtime feedback scan
% const MdhBitField MDH_HPFEEDBACK        (2UL);       ///< High perfomance feedback scan
% const MdhBitField MDH_ONLINE            (3UL);       ///< processing should be done online
% const MdhBitField MDH_OFFLINE           (4UL);       ///< processing should be done offline
% const MdhBitField MDH_SYNCDATA          (5UL);       ///< readout contains synchroneous data
% const MdhBitField MDH_LASTSCANINCONCAT  (8UL);       ///< Flag for last scan in concatination
% 
% const MdhBitField MDH_RAWDATACORRECTION (10UL);      ///< Correct the rawadata with the rawdata correction factor
% const MdhBitField MDH_LASTSCANINMEAS    (11UL);      ///< Flag for last scan in measurement
% const MdhBitField MDH_SCANSCALEFACTOR   (12UL);      ///< Flag for scan specific additional scale factor
% const MdhBitField MDH_2NDHADAMARPULSE   (13UL);      ///< 2nd RF exitation of HADAMAR
% const MdhBitField MDH_REFPHASESTABSCAN  (14UL);      ///< reference phase stabilization scan
% const MdhBitField MDH_PHASESTABSCAN     (15UL);      ///< phase stabilization scan
% const MdhBitField MDH_D3FFT             (16UL);      ///< execute 3D FFT
% const MdhBitField MDH_SIGNREV           (17UL);      ///< sign reversal
% const MdhBitField MDH_PHASEFFT          (18UL);      ///< execute phase fft
% const MdhBitField MDH_SWAPPED           (19UL);      ///< swapped phase/readout direction
% const MdhBitField MDH_POSTSHAREDLINE    (20UL);      ///< shared line
% const MdhBitField MDH_PHASCOR           (21UL);      ///< phase correction data
% const MdhBitField MDH_PATREFSCAN        (22UL);      ///< additonal scan for PAT reference line/partition
% const MdhBitField MDH_PATREFANDIMASCAN  (23UL);      ///< additonal scan for PAT reference line/partition that is also used as image scan
% const MdhBitField MDH_REFLECT           (24UL);      ///< reflect line
% const MdhBitField MDH_NOISEADJSCAN      (25UL);      ///< noise adjust scan 
% const MdhBitField MDH_SHARENOW          (26UL);      ///< all lines are acquired from the actual and previous e.g. phases
% const MdhBitField MDH_LASTMEASUREDLINE  (27UL);      ///< indicates that the current line is the last measured line of all succeeding e.g. phases
% const MdhBitField MDH_FIRSTSCANINSLICE  (28UL);      ///< indicates first scan in slice (needed for time stamps)
% const MdhBitField MDH_LASTSCANINSLICE   (29UL);      ///< indicates  last scan in slice (needed for time stamps)
% const MdhBitField MDH_TREFFECTIVEBEGIN  (30UL);      ///< indicates the begin time stamp for TReff (triggered measurement)
% const MdhBitField MDH_TREFFECTIVEEND    (31UL);      ///< indicates the   end time stamp for TReff (triggered measurement)
% const MdhBitField MDH_MDS_REF_POSITION  (32UL);      ///< indicates the reference position for move during scan images (must be set once per slice/partition in MDS mode)
% const MdhBitField MDH_SLC_AVERAGED      (33UL);      ///< indicates avveraged slice for slice partial averaging scheme
% const MdhBitField MDH_TAGFLAG1          (34UL);      ///< adjust scan 
% 
% const MdhBitField MDH_CT_NORMALIZE              (35UL);  ///< Marks scans used to calculate correction maps for TimCT-Prescan normalize
% const MdhBitField MDH_SCAN_FIRST                (36UL);  ///< Marks the first scan of a particular map
% const MdhBitField MDH_SCAN_LAST                 (37UL);  ///< Marks the last scan of a particular map
% 
% const MdhBitField MDH_FIRST_SCAN_IN_BLADE       (40UL);  ///< Marks the first line of a blade
% const MdhBitField MDH_LAST_SCAN_IN_BLADE        (41UL);  ///< Marks the last line of a blade
% const MdhBitField MDH_LAST_BLADE_IN_TR          (42UL);  ///< Set for all lines of the last BLADE in each TR interval
%                                                           
% const MdhBitField MDH_PACE                      (44UL);  ///< Distinguishes PACE scans from non PACE scans.
%                                                           
% const MdhBitField MDH_RETRO_LASTPHASE           (45UL);  ///< Marks the last phase in a heartbeat
% const MdhBitField MDH_RETRO_ENDOFMEAS           (46UL);  ///< Marks an ADC at the end of the measurement
% const MdhBitField MDH_RETRO_REPEATTHISHEARTBEAT (47UL);  ///< Repeat the current heartbeat when this bit is found
% const MdhBitField MDH_RETRO_REPEATPREVHEARTBEAT (48UL);  ///< Repeat the previous heartbeat when this bit is found
% const MdhBitField MDH_RETRO_ABORTSCANNOW        (49UL);  ///< Just abort everything
% const MdhBitField MDH_RETRO_LASTHEARTBEAT       (50UL);  ///< This adc is from the last heartbeat (a dummy)
% const MdhBitField MDH_RETRO_DUMMYSCAN           (51UL);  ///< This adc is just a dummy scan, throw it away
% const MdhBitField MDH_RETRO_ARRDETDISABLED      (52UL);  ///< Disable all arrhythmia detection when this bit is found


function out = PMU_training(PMU);

tmp = []; tmp_time = [];
for i1 = 1:size(PMU.Data,2);
    if ~isempty(findstr(PMU.PackedID{i1},'PMULearnPhase'));
        tmp = [tmp,PMU.Data{1,i1}'];
        len = length(PMU.Data{1,i1});
        pmu_timestamp = PMU.TimeStamp(i1,2); % lower 32 bit word time stamp for j-th packet
        tmp_time = [tmp_time, pmu_timestamp*2.5e-3+ 2.5e-3*[0:len-1]]; % ECG samples are in 0.1ms increments
        % from pmu timestamp which is in 2.5e-3 increments from 00:00:00
    end;
end;
ecg1 = bitand(tmp,32767);% lower 16 bits
trigger1 = bitand(tmp, 2^30);
mask1 = bitand(tmp, hex2dec('FFFF0000'))/2^16; % upper 16 bits
ecg1_time = tmp_time;

tmp = [];for i1 = 1:size(PMU.Data,2); if ~isempty(findstr(PMU.PackedID{i1},'PMULearnPhase'));tmp = [tmp,PMU.Data{2,i1}'];end;end;
ecg2 = bitand(tmp,32767);
trigger2 = bitand(tmp, 2^30);

tmp = [];for i1 = 1:size(PMU.Data,2); if ~isempty(findstr(PMU.PackedID{i1},'PMULearnPhase'));tmp = [tmp,PMU.Data{3,i1}'];end;end;
ecg3 = bitand(tmp,32767);
trigger3 = bitand(tmp, 2^30);

tmp = [];for i1 = 1:size(PMU.Data,2); if ~isempty(findstr(PMU.PackedID{i1},'PMULearnPhase'));tmp = [tmp,PMU.Data{4,i1}'];end;end;
ecg4 = bitand(tmp,32767);
trigger4 = bitand(tmp, 2^30);

tmp = [];for i1 = 1:size(PMU.Data,2); if ~isempty(findstr(PMU.PackedID{i1},'PMULearnPhase'));tmp = [tmp,PMU.Data{5,i1}'];end;end;
puls = bitand(tmp,32767);
trigger_puls = bitand(tmp, 2^29);

tmp = [];for i1 = 1:size(PMU.Data,2); if ~isempty(findstr(PMU.PackedID{i1},'PMULearnPhase'));tmp = [tmp,PMU.Data{6,i1}'];end;end;
resp = bitand(tmp,32767);
trigger_resp = bitand(tmp, 2^28);


ecg1 = (ecg1-2048); 
ecg2 = (ecg2-2048);
ecg3 = (ecg3-2048);
ecg4 = (ecg4-2048);

out.ecg1=ecg1;
out.trigger_ecg1 = trigger1;
out.mask1 = mask1;
out.ecg1_time = ecg1_time;
out.ecg2=ecg2;
out.trigger_ecg2 = trigger2;
out.ecg3=ecg3;
out.trigger_ecg3 = trigger3;
out.ecg4=ecg4;
out.trigger_ecg4 = trigger4;
out.puls=puls;
out.trigger_puls = trigger_puls;
out.puls=resp;
out.trigger_resp = trigger_resp;

