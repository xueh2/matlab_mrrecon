% Try to read the raw data from the scanner
% function [Data, asc, prot] = Read_RawData(fname)

function [Data, asc, prot] = Read_RawData(fname)
%function [Data, asc, prot, header] = Read_RawData(fname)
% asc(i).ushUsedChannels: Total # of channels, fixed
% asc(i).ushChannelId: Channel #

tic
fid = fopen(fname);
L_P = fread(fid, 1, 'uint32');
prot = fread(fid, L_P-4, 'uint8');
% c0 = clock;
fseek(fid, 7*4,'cof');
Samples = fread(fid, 1, 'uint16');
STATUS = fseek(fid, 128-30,'cof');
STATUS = fseek(fid, 8*Samples + 4, 'cof');
i = 1;
maxCol = Samples;
while (STATUS+1)
    i = i + 1;
    fseek(fid, 7*4 - 4 ,'cof');
    Samples = fread(fid, 1, 'uint16');
    STATUS = fseek(fid, 128-30,'cof');
    STATUS = fseek(fid, 8*Samples + 4, 'cof'); % + 1 to test if it is end of the file
    if ( Samples > maxCol )
        maxCol = Samples;
    end
end
% The above is to detect the number of lines
N = i,

fseek(fid,0,-1) % "rewinds" the file.
STATUS = fseek(fid, L_P, 'cof');
fseek(fid, 7*4,'cof');
Samples = fread(fid, 1, 'uint16'),
STATUS = fseek(fid, -30,'cof');

asc0 = struct( ...
'ulFlagsAndDMALength',  0, ...
'lMeasUID', 0, ...
'ulScanCounter',  0,...
'ulTimeStamp',  0,...
'ulPMUTimeStamp',  0,...
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
'aushIceProgramPara',  zeros(4,1),...
'aushFreePara',  zeros(4,1),...
'sSD',  zeros(14,1),...
'ushChannelId',  0,...
'ushPTABPosNeg',  0 ...
);

asc = repmat(asc0,1,N);
Data = complex(single(zeros(maxCol,N)), single(zeros(maxCol,N)));
k = 0;
for i=1:N
       
    if ( mod(i, 1e4) == 0 )
        disp([num2str(i) ' lines read ... ']);
    end
    
    data = fread(fid, 128,  'uint8');
    data = uint8(data);
    lenData = numel(data);
    ind = 1;
    asc(i).ulFlagsAndDMALength = typecast(data(1:4),            'uint32'); ind = ind + 4; if ( ind >= lenData ) break; end
    asc(i).lMeasUID = typecast(data(ind:ind+3),                 'int32');  ind = ind + 4; if ( ind >= lenData ) break; end
    asc(i).ulScanCounter = typecast(data(ind:ind+3),            'uint32'); ind = ind + 4; if ( ind >= lenData ) break; end
    asc(i).ulTimeStamp = typecast(data(ind:ind+3),              'uint32'); ind = ind + 4; if ( ind >= lenData ) break; end
    asc(i).ulPMUTimeStamp = typecast(data(ind:ind+3),           'uint32'); ind = ind + 4; if ( ind >= lenData ) break; end
    asc(i).aulEvalInfoMask = typecast(data(ind:ind+7),          'uint32'); ind = ind + 8; if ( ind >= lenData ) break; end
    asc(i).ushSamplesInScan = typecast(data(ind:ind+1),         'uint16'); ind = ind + 2; if ( ind >= lenData ) break; end
    asc(i).ushUsedChannels = typecast(data(ind:ind+1),          'uint16'); ind = ind + 2; if ( ind >= lenData ) break; end
    asc(i).sLC = typecast(data(ind:ind+27),                     'uint16'); ind = ind + 2*14; if ( ind >= lenData ) break; end
    asc(i).sCutOff = typecast(data(ind:ind+3),                  'single');  ind = ind + 4; if ( ind >= lenData ) break; end
    asc(i).ushKSpaceCentreColumn = typecast(data(ind:ind+1),    'uint16'); ind = ind + 2; if ( ind >= lenData ) break; end
    asc(i).ushCoilSelect = typecast(data(ind:ind+1),            'uint16'); ind = ind + 2; if ( ind >= lenData ) break; end
    asc(i).fReadOutOffcentre = typecast(data(ind:ind+3),        'single'); ind = ind + 4; if ( ind >= lenData ) break; end
    asc(i).ulTimeSinceLastRF = typecast(data(ind:ind+3),        'uint32'); ind = ind + 4; if ( ind >= lenData ) break; end
    asc(i).ushKSpaceCentreLineNo = typecast(data(ind:ind+1),    'uint16'); ind = ind + 2; if ( ind >= lenData ) break; end
    asc(i).ushKSpaceCentrePartitionNo = typecast(data(ind:ind+1), 'uint16'); ind = ind + 2; if ( ind >= lenData ) break; end
    asc(i).aushIceProgramPara = typecast(data(ind:ind+7),       'uint16'); ind = ind + 2*4; if ( ind >= lenData ) break; end
    asc(i).aushFreePara = typecast(data(ind:ind+7),             'uint16'); ind = ind + 2*4; if ( ind >= lenData ) break; end
    asc(i).sSD = typecast(data(ind:ind+27),                     'uint16'); ind = ind + 2*14; if ( ind >= lenData ) break; end
    asc(i).ushChannelId = typecast(data(ind:ind+1),             'uint16'); ind = ind + 2; if ( ind >= lenData ) break; end
    asc(i).ushPTABPosNeg = typecast(data(ind:ind+1),            'uint16'); ind = ind + 2; 
    
%     fseek(fid, -128,'cof')
%     
%     asc(i).ulFlagsAndDMALength =        fread(fid, 1, 'uint32');
%     asc(i).lMeasUID =                   fread(fid, 1, 'int32');
%     asc(i).ulScanCounter =              fread(fid, 1, 'uint32');
%     asc(i).ulTimeStamp  =               fread(fid, 1, 'uint32');
%     asc(i).ulPMUTimeStamp =             fread(fid, 1, 'uint32');
%     asc(i).aulEvalInfoMask =            fread(fid, 2, 'uint32');
%     asc(i).ushSamplesInScan =           fread(fid, 1, 'uint16');
%     asc(i).ushUsedChannels =            fread(fid, 1, 'uint16'); % Channels % + 8 =8
%     asc(i).sLC =                        fread(fid, 14, 'uint16')'; % + 7 = 15
%     asc(i).sCutOff =                    fread(fid, 1, 'float'); % + 1 = 17
%     asc(i).ushKSpaceCentreColumn =      fread(fid, 1, 'uint16');
%     asc(i).ushCoilSelect =              fread(fid, 1, 'uint16');
%     asc(i).fReadOutOffcentre =          fread(fid, 1, 'float');
%     asc(i).ulTimeSinceLastRF =          fread(fid, 1, 'uint32'); % +3 = 19
%     asc(i).ushKSpaceCentreLineNo =      fread(fid, 1, 'uint16');
%     asc(i).ushKSpaceCentrePartitionNo = fread(fid, 1, 'uint16');
%     asc(i).aushIceProgramPara =         fread(fid, 4, 'uint16'); % + 3 = 22
%     asc(i).aushFreePara =               fread(fid, 4, 'uint16') ; % +2
%     asc(i).sSD =                        fread(fid, 14, 'uint16'); % + 7 = 31, No. 5 is No. of lines
%     asc(i).ushChannelId =               fread(fid, 1, 'uint16');
%     asc(i).ushPTABPosNeg =              fread(fid, 1, 'uint16'); % +1 = 32;
    t = 0; [t, count] =                 fread(fid, 2*asc(i).ushSamplesInScan, 'float32');
    Data(1:asc(i).ushSamplesInScan,i) = single(complex( t(1:2:end), t(2:2:end) ));
end
fclose(fid)
toc