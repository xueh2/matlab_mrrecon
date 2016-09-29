% Try to read the write the raw data from the scanner
% function [Data, asc, prot] = Read_RawData(fname)
 
function measData = Read_RawData_VD11A(fname)
% fname: file name of the meas.dat
% flagIsPostVD11A: is the file in VD11A (or beyond) format?

%function [Data, asc, prot, header] = Read_RawData(fname)
% asc(i).ushUsedChannels: Total # of channels, fixed
% asc(i).ushChannelId: Channel #
 
fid = fopen(fname);

hdSize = fread(fid, 1, 'uint32'); % doesn't seem to be important
measNum = fread(fid, 1, 'uint32'); % # of measurements in this file
% start reading MrParcRaidFileEntry (maximum is 64)
measCount = 0;
measOffset = zeros(measNum, 1);
measLength = zeros(measNum, 1);
while measCount < measNum
    measId = fread(fid, 1, 'uint32');
    if measId > 0 % there is a measurement (adjustments, pre-scan, noise-scan, "real" scan...)
        measCount = measCount + 1;
        fileId = fread(fid, 1, 'uint32'); % no use as of 3/22/2011
        measOffset(measCount) = fread(fid, 1, 'uint64'); % offset from the bof
        measLength(measCount) = fread(fid, 1, 'uint64');
        patientName = char(fread(fid, 64, 'char')); % no use 3/22/2011
        protocolName = char(fread(fid, 64, 'char')); % no use 3/22/2011
    end
end
% now the MrParcRaidFileEntry is parsed, use cell structure to hold all 
% measurements
measData = cell(measNum, 1);
N = zeros(measNum, 1);
totalLine = zeros(measNum, 1);
for mi= 1 : measNum
    fseek(fid, measOffset(mi), 'bof');
    
    % start reading header
    L_P = fread(fid, 1, 'uint32');
    measData{mi}.prot = fread(fid, L_P-4, 'uint8');
    i = 0;
    totalLineCount = 0;
    notAlignBytes = true; %fseek(fid, 4,'cof');% seek the first 4 bytes for status
    while (notAlignBytes)
        % The scan header (first 4 bytes has been read), only one such header for all channels
        fseek(fid, 10*4,'cof');
        evalInfoMask = fread(fid, 2, 'uint32');
        Samples = fread(fid, 1, 'uint16');
        channelNum = fread(fid, 1, 'uint16');
        fseek(fid, 140,'cof'); % 140 = 192-52
        for ci = 1 : channelNum
                fseek(fid, 32,'cof'); % ch header
                fseek(fid, 8*Samples, 'cof');
                totalLineCount = totalLineCount + 1;
        end
        i = i + 1;
        if evalInfoMask(1) == 1 % AcqEnd
            notAlignBytes = false;
        end
        
    end
    % The above is to detect the number of lines per channel
    N(mi) = i;
    totalLine(mi) = totalLineCount;
end

 
asc0 = struct( ...
'ulFlagsAndDMALength',  0, ...
'lMeasUID', 0, ...
'ulScanCounter',  0,...
'ulTimeStamp',  0,...
'ulPMUTimeStamp',  0,...
'SystemType', 0, ...
'PTABPosDelay', 0, ...
'PTABPosX', 0, ...
'PTABPosY', 0, ...
'PTABPosZ', 0, ...
'Reserved1', 0, ...
'aulEvalInfoMask',  zeros(2,1),...
'ushSamplesInScan',  0,...
'ushUsedChannels',  0,...
'sLC',  zeros(14,1),...
'sCutOff',  0,...
'ushKSpaceCentreColumn',  0,...
'ushCoilSelect',  0,...
'fReadOutOffcentre',  0,...
'ulTimeSinceLastRF',  0,...
'ushKSpaceCentreLineNo',  0,...
'ushKSpaceCentrePartitionNo',  0,...
'sSD',  zeros(14,1),...
'aushIceProgramPara',  zeros(4,1),...
'ReservedPara', zeros(4, 1), ...
'ApplicationCounter', 0, ...
'ApplicationMask', 0, ...
'CRC', 0, ...
'ushChannelId', -ones(128,1, 'int16') ...
);
 
for mi= 1 : measNum

    totalLineCount = 1;
    asc = repmat(asc0,1,N(mi));
    Data = complex(single(zeros(Samples,totalLine(mi))), single(zeros(Samples,totalLine(mi))));
    
    fseek(fid, measOffset(mi), 'bof');
    L_P = fread(fid, 1, 'uint32');
    fseek(fid, L_P-4, 'cof');

    for i=1:N(mi)
        %i = i + 1; % if mod(i,1000)==0, i=i,clock, end
        asc(i).ulFlagsAndDMALength =        fread(fid, 1, 'uint32');
        asc(i).lMeasUID =                   fread(fid, 1, 'int32');
        asc(i).ulScanCounter =              fread(fid, 1, 'uint32');
        asc(i).ulTimeStamp  =               fread(fid, 1, 'uint32');
        asc(i).ulPMUTimeStamp =             fread(fid, 1, 'uint32');
        asc(i).SystemType =                 fread(fid, 1, 'uint16');
        asc(i).PTABPosDelay =               fread(fid, 1, 'uint16');
        asc(i).PTABPosX =                   fread(fid, 1, 'single');
        asc(i).PTABPosY =                   fread(fid, 1, 'single');
        asc(i).PTABPosZ =                   fread(fid, 1, 'single');
        asc(i).Reserved1 =                  fread(fid, 1, 'uint32');
        asc(i).aulEvalInfoMask =            fread(fid, 2, 'uint32');
        asc(i).ushSamplesInScan =           fread(fid, 1, 'uint16');
        asc(i).ushUsedChannels =            fread(fid, 1, 'uint16'); % Channels % + 8 =8
        asc(i).sLC =                        fread(fid, 14, 'uint16')'; % + 7 = 15
        asc(i).sCutOff =                    fread(fid, 1, 'single'); % + 1 = 17
        asc(i).ushKSpaceCentreColumn =      fread(fid, 1, 'uint16');
        asc(i).ushCoilSelect =              fread(fid, 1, 'uint16');
        asc(i).fReadOutOffcentre =          fread(fid, 1, 'single');
        asc(i).ulTimeSinceLastRF =          fread(fid, 1, 'uint32'); % +3 = 19
        asc(i).ushKSpaceCentreLineNo =      fread(fid, 1, 'uint16');
        asc(i).ushKSpaceCentrePartitionNo = fread(fid, 1, 'uint16');
        asc(i).sSD =                        fread(fid, 14, 'uint16');
        asc(i).aushIceProgramPara =         fread(fid, 24, 'uint16');
        asc(i).ReservedPara =               fread(fid, 4, 'uint16') ;
        asc(i).ApplicationCounter =         fread(fid, 1, 'uint16');
        asc(i).ApplicationMask =            fread(fid, 1, 'uint16');
        asc(i).CRC =                        fread(fid, 1, 'uint32');
        
%         asc(i).ushChannelId = zeros(asc(i).ushUsedChannels, 1);
        % Only extract ChannelId from the channel header

        for ci= 1 : asc(i).ushUsedChannels
            fseek(fid, 4*6, 'cof');
            asc(i).ushChannelId(ci) =       fread(fid, 1, 'uint16');
            fseek(fid, 2+4, 'cof');
            t = 0; [t, count] =             fread(fid, 2*asc(i).ushSamplesInScan, 'float32');
            temp = asc(i).sLC; %Slice_Counter  = temp(3);
            Data(1:asc(i).ushSamplesInScan,totalLineCount) = single(complex( t(1:2:end), t(2:2:end) ));
            totalLineCount = totalLineCount + 1;
        end

    end
    measData{mi}.asc = asc;
    measData{mi}.Data = Data;
end
fclose(fid);