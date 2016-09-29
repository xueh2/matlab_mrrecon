function [headers, Data, asc, dataH5] = readSiemensDataVD(h5filename, scanno)
% read in the siemens meas dat file header for VD

% the header
groupName = '/files/';
groupName = [groupName num2str(scanno) '/MeasurementHeader'];

h5headers = h5read(h5filename, groupName);
    
headers.Config = char(h5headers.buffers{1}.buf_{1});
headers.Config = permute(headers.Config, [2 1]);

headers.Dicom = char(h5headers.buffers{1}.buf_{2});
headers.Dicom = permute(headers.Dicom, [2 1]);

headers.MrcTime = char(h5headers.buffers{1}.buf_{3});
headers.MrcTime = permute(headers.MrcTime, [2 1]);

headers.MeasYaps = char(h5headers.buffers{1}.buf_{4});
headers.MeasYaps = permute(headers.MeasYaps, [2 1]);

headers.Phoenix = char(h5headers.buffers{1}.buf_{5});
headers.Phoenix = permute(headers.Phoenix, [2 1]);

headers.Spice = char(h5headers.buffers{1}.buf_{6});
headers.Spice = permute(headers.Spice, [2 1]);

[path, name, ext] = fileparts(h5filename);

fid = fopen(fullfile(path, [name '_Config.txt']), 'w');
fprintf(fid, '%s', char(headers.Config));
fclose(fid);

fid = fopen(fullfile(path, [name '_Dicom.txt']), 'w');
fprintf(fid, '%s', char(headers.Dicom));
fclose(fid);

fid = fopen(fullfile(path, [name '_MeasYaps.txt']), 'w');
fprintf(fid, '%s', char(headers.MeasYaps));
fclose(fid);

fid = fopen(fullfile(path, [name '_Phoenix.txt']), 'w');
fprintf(fid, '%s', char(headers.Phoenix));
fclose(fid);

fid = fopen(fullfile(path, [name '_Spice.txt']), 'w');
fprintf(fid, '%s', char(headers.Spice));
fclose(fid);

% the raw data

groupName = '/files/';
groupName = [groupName num2str(scanno) '/data'];

data = h5read(h5filename, groupName);
dataH5 = data;

N = numel(data.data) - 1; % the last scan is the stop scan

if ( data.scanHeader.aulEvalInfoMask(1, end) ~= 1 )
    warning('Incomplete scan found ... ');
    N = N + 1;
end

asc0 = struct( ...
'ulFlagsAndDMALength',  0, ...
'lMeasUID', 0, ...
'ulScanCounter',  0,...
'ulTimeStamp',  0,...
'ulPMUTimeStamp',  0,...
'ushSystemType', 0,...
'ulPTABPosDelay', 0,...
'lPTABPosX', 0,...
'lPTABPosY', 0,...
'lPTABPosZ', 0,...
'ulReserved1', 0,...
'aulEvalInfoMask',  zeros(2,1),...
'ushSamplesInScan',  0,...
'ushUsedChannels',  0,...
'sLC',  zeros(14,1),...
'sCutOff',  0,...
'ushChannelId',  0,...
'ushKSpaceCentreColumn',  0,...
'ushCoilSelect',  0,...
'fReadOutOffcentre',  0,...
'ulTimeSinceLastRF',  0,...
'ushKSpaceCentreLineNo',  0,...
'ushKSpaceCentrePartitionNo',  0,...
'sSliceData_sSlicePosVec', 0, ...
'sSliceData_aflQuaternion', 0, ...
'aushIceProgramPara',  zeros(4,1),...
'aushReservedPara', 0,...
'ushApplicationCounter', 0,...
'ushApplicationMask', 0,...
'ulCRC', 0 );

% find the maxCol
maxCol = max(data.scanHeader.ushSamplesInScan(:));
nCha = numel(data.data{1}.data);

asc = repmat(asc0,1,N*nCha);
% Data = complex(single(zeros(maxCol,N*nCha)), single(zeros(maxCol,N*nCha)));
Data = complex(single(zeros(maxCol,nCha,N)), single(zeros(maxCol,nCha,N)));

for scan=1:N
    disp(['Read scan : ' num2str(scan) ' ... ']);
    for cha=1:nCha
        
        ii = (scan-1)*nCha+cha;
        sampleLen = data.scanHeader.ushSamplesInScan(scan);

        if ( cha==1 )
            asc(ii).ulFlagsAndDMALength = data.scanHeader.ulFlagsAndDMALength(scan);
            asc(ii).lMeasUID = data.scanHeader.lMeasUID(scan);
            asc(ii).ulScanCounter = data.scanHeader.ulScanCounter(scan);
            asc(ii).ulTimeStamp = data.scanHeader.ulTimeStamp(scan);
            asc(ii).ulPMUTimeStamp = data.scanHeader.ulPMUTimeStamp(scan);
            asc(ii).ushSystemType = data.scanHeader.ushSystemType(scan);
            asc(ii).ulPTABPosDelay = data.scanHeader.ulPTABPosDelay(scan);
            asc(ii).lPTABPosX = data.scanHeader.lPTABPosX(scan);
            asc(ii).lPTABPosY = data.scanHeader.lPTABPosY(scan);
            asc(ii).lPTABPosZ = data.scanHeader.lPTABPosZ(scan);
            asc(ii).ulReserved1 = data.scanHeader.ulReserved1(scan);

            asc(ii).aulEvalInfoMask = data.scanHeader.aulEvalInfoMask(:,scan);
            asc(ii).ushSamplesInScan = data.scanHeader.ushSamplesInScan(scan);
            asc(ii).ushUsedChannels = data.scanHeader.ushUsedChannels(scan);

            asc(ii).sLC = [data.scanHeader.sLC.ushLine(scan) ...
                           data.scanHeader.sLC.ushAcquisition(scan) ...
                           data.scanHeader.sLC.ushSlice(scan) ...
                           data.scanHeader.sLC.ushPartition(scan) ...
                           data.scanHeader.sLC.ushEcho(scan) ...
                           data.scanHeader.sLC.ushPhase(scan) ...
                           data.scanHeader.sLC.ushRepetition(scan) ...
                           data.scanHeader.sLC.ushSet(scan) ...
                           data.scanHeader.sLC.ushSeg(scan) ...
                           data.scanHeader.sLC.ushIda(scan) ...
                           data.scanHeader.sLC.ushIdb(scan) ...
                           data.scanHeader.sLC.ushIdc(scan) ...
                           data.scanHeader.sLC.ushIdd(scan) ...
                           data.scanHeader.sLC.ushIde(scan) ];

            asc(ii).ushChannelId = cha - 1;

            asc(ii).sCutOff = [data.scanHeader.sCutOff.ushPre(scan) data.scanHeader.sCutOff.ushPost(scan)];

            asc(ii).ushKSpaceCentreColumn = data.scanHeader.ushKSpaceCentreColumn(scan);
            asc(ii).ushCoilSelect = data.scanHeader.ushCoilSelect(scan);
            asc(ii).fReadOutOffcentre = data.scanHeader.fReadOutOffcentre(scan);
            asc(ii).ulTimeSinceLastRF = data.scanHeader.ulTimeSinceLastRF(scan);
            asc(ii).ushKSpaceCentreLineNo = data.scanHeader.ushKSpaceCentreLineNo(scan);
            asc(ii).ushKSpaceCentrePartitionNo = data.scanHeader.ushKSpaceCentrePartitionNo(scan);

            asc(ii).sSliceData_sSlicePosVec = [data.scanHeader.sSliceData.sSlicePosVec.flCor(scan) data.scanHeader.sSliceData.sSlicePosVec.flSag(scan) data.scanHeader.sSliceData.sSlicePosVec.flTra(scan)];
            asc(ii).sSliceData_aflQuaternion = data.scanHeader.sSliceData.aflQuaternion(:, scan);

            asc(ii).aushIceProgramPara = data.scanHeader.aushIceProgramPara(:,scan);
            asc(ii).aushReservedPara = data.scanHeader.aushReservedPara(:,scan);
            asc(ii).ushApplicationCounter = data.scanHeader.ushApplicationCounter(scan);
            asc(ii).ushApplicationMask = data.scanHeader.ushApplicationMask(scan);
            asc(ii).ulCRC = data.scanHeader.ulCRC(scan);
        else
            asc(ii) = asc((scan-1)*nCha+1);
            asc(ii).ushChannelId = cha - 1;
        end
        
        try
            cData = complex(data.data{scan}.data{cha}.real, data.data{scan}.data{cha}.imag);
        catch
            cData = zeros(sampleLen, 1);
        end
        
        Data(1:sampleLen, cha, scan) = cData;
    end
end
