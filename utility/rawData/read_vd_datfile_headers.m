function [MrParcRaidFileHeader, MrParcRaidFileEntry, headers, protocol]=read_vd_datfile_headers(datfilename);
% function [MrParcRaidFileHeader, MrParcRaidFileEntry, headers, protocol]=read_vd_datfile_headers(datfilename);
%
% function to read Siemens dat file headers for meas.dat (for MultiRaid
% files found in VD11 and higher)
%
% headers{i} is output structure for multiraid dataset i containing ascii protocol data in raw text form
%     headers{i}.Config
%     headers{i}.Dicom
%     headers{i}.MRCTime
%     headers{i}.Meas
%     headers{i}.MeasYaps
%     headers{i}.Phoenix
%     headers{i}.Spice

% protocol (optional output) is Matlab structure corresponding to MeasYaps
       
%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  NIH NHLBI                          *
%     ***************************************

fid=fopen(datfilename,'r');

MrParcRaidFileHeader.hdsize=fread(fid,1,'uint32');
switch MrParcRaidFileHeader.hdsize
    % MR_PARC_RAID_ALLDATA: normal measurements
    % MR_PARC_RAID_MDHONLY: measurements without raw data
    % MR_PARC_RAID_HDONLY : measurements without mdh and raw data
    % MR_PARC_RAID_LEGACY_THR: pre VD11A raidfile
    case 0
        MrParcRaidFileHeader.hdsize = 'MR_PARC_RAID_ALLDATA';
    case 1
        MrParcRaidFileHeader.hdsize = 'MR_PARC_RAID_MDHONLY';
    case 2
        MrParcRaidFileHeader.hdsize = 'MR_PARC_RAID_HDONLY';
end
MrParcRaidFileHeader.count=fread(fid,1,'uint32');
for count = 1:MrParcRaidFileHeader.count
    MrParcRaidFileEntry(count).measID=fread(fid,1,'uint32');
    MrParcRaidFileEntry(count).fileID=fread(fid,1,'uint32');
    MrParcRaidFileEntry(count).off=fread(fid,1,'uint64');
    MrParcRaidFileEntry(count).len=fread(fid,1,'uint64');
    MrParcRaidFileEntry(count).patName=char(fread(fid,64,'char'))';
    MrParcRaidFileEntry(count).protName=char(fread(fid,64,'char'))';
end

if nargout < 3; return; end

for count = 1:MrParcRaidFileHeader.count
    fseek(fid,MrParcRaidFileEntry(count).off,-1);% skip to file offset from bof
    dma_len=fread(fid,1,'int32');
    nbuffers=fread(fid,1,'int32');
    for i=1:nbuffers
       tmp=1;
       name{i}=[];
       while tmp~=0;
           tmp=fread(fid,1,'char');
           name{i}=[name{i},tmp];
       end
       name{i}=char(name{i}(1:end-1));
       length(i)=fread(fid,1,'int32');
       tmp=char(fread(fid,length(i),'char')');
       eval(['headers{count}.',name{i},'= tmp;']);   
    end

end

if nargout==4
    for count = 1:MrParcRaidFileHeader.count
        protocol{count}=MeasYaps2struct(headers{count});
    end
end

fclose(fid);

return


function header=MeasYaps2struct(headers);
% function header=MeasYaps2struct(headers);
%
% function to covert MeasYaps from text to Matlab structure

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  NIH NHLBI                          *
%     ***************************************

MeasYaps=headers.MeasYaps;
lf = find(double(MeasYaps)==10); % find positions of line feeds
lf = [0,lf];

for line=1:length(lf)-1
    tline = MeasYaps(lf(line)+1:lf(line+1)-1);
    if ~ischar(tline); return; end
    index=findstr(tline,'=');
    ASCCONV =findstr(tline,'ASCCONV');
    if ~isempty(index) & isempty(ASCCONV)
        fieldname=deblank(tline(1:findstr(tline,'=')-1));
        fieldvalue=deblank(tline(findstr(tline,'=')+1:end));
        if max(isletter(fieldvalue))==0
            tmpfieldvalue=str2num(fieldvalue);
            if ~isempty(tmpfieldvalue);
                fieldvalue = tmpfieldvalue;
            end
        end
        fieldname=strrep(fieldname,'[','{');
        fieldname=strrep(fieldname,']','}');
        index2=findstr(fieldname,'{');
        index3=findstr(fieldname,'}');


        eb = size(index2,2);
        ee = size(index3,2);
      
        strs = eb + 1;
        nums = eb;
        
        if ~isempty(index2)

            %split out the character portions of the statement
            str(1) = {fieldname(1:index2(1))};
            for i = 2:strs-1
                str(i) = {fieldname(index3(i-1):index2(i))};
            end
            str(strs) = {fieldname(index3(strs-1):end)};
            
			%pull out, increment and reconvert to string the numeric portion of the
			%statement
			for i=1:nums
                numstr(i) = {num2str(str2num(fieldname(index2(i)+1:index3(i)-1))+1)};
			end

			%combine the string and the numeric portions of the statement
			nstr = '';
            for i=1:eb
                nstr = strcat(nstr, str(i), numstr(i));
			end
            
           %convert from cell array to char array before pushing back in.
           fieldname = char(strcat(nstr,str(strs)));
        end
        fieldname = strrep(fieldname,'_','');
        curlybracket = findstr(fieldname,'{');
        if ~isempty(curlybracket)
            for b = 1:length(curlybracket)
               if ~iscell(['header.',fieldname(1:curlybracket(b)-1)])
                   fieldname = [fieldname(1:curlybracket(b)-1),'_',fieldname(curlybracket(b):end)];
                   curlybracket = findstr(fieldname,'{');
               end
            end
        end
        eval(['header.',fieldname,'= fieldvalue;']);
    end
end

function header=DatHeader2struct(datheader);
% function header=DatHeader2struct(datheader);
%
% function to convert DatHeader from text to Matlab structure

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  NIH NHLBI                          *
%     ***************************************

lf = find(double(datheader)==10); % find positions of line feeds
lf = [0,lf];

for line=1:length(lf)-1
    tline = datheader(lf(line)+1:lf(line+1)-1);
    if ~ischar(tline); return; end
    index=findstr(tline,'ParamString."');    
    if ~isempty(index)
        % example        <ParamString."tSequenceVariant">  { "SK\SP"  }
        index_quote=findstr(tline,'"');
        if length(index_quote)>=4
            fieldname=deblank(tline(index_quote(1)+1:index_quote(2)-1));
            fieldvalue=deblank(tline(index_quote(3)+1:index_quote(4)-1));
            eval(['header.',fieldname,'= fieldvalue;']);
        end
    end
    index=findstr(tline,'ParamLong."');
    if ~isempty(index)
        index_quote=findstr(tline,'"');
        if length(index_quote)>=2
            fieldname=deblank(tline(index_quote(1)+1:index_quote(2)-1));
        end
        if ~isempty(fieldname)
            % case 1....value is contain between { } on same line
            %    for example:     <ParamLong."MeasUID">  { 71  }
            start = findstr(tline,'{')+1;
            stop  = findstr(tline,'}')-1;
            fieldvalue=str2num(deblank(tline(start:stop)));
            eval(['header.',fieldname,'= fieldvalue;']);        
             % case 2....value is between {} on next lines
             %   for example
                %       <ParamLong."NoOfFourierLines"> 
                %       {
                %         120 
                %       }
            if isempty(start)
                line=line+1;
                tline = datheader(lf(line)+1:lf(line+1)-1);
                if ~isempty (findstr(tline,'{'))
                    line=line+1;
                    tline = datheader(lf(line)+1:lf(line+1)-1);
                    if isempty(findstr(tline,'<Comment>')) % skip comment line
                        fieldvalue=str2num(deblank(tline));
                        eval(['header.',fieldname,'= fieldvalue;']);
                    else
                        line=line+1;
                        tline = datheader(lf(line)+1:lf(line+1)-1);
                        fieldvalue=str2num(deblank(tline));
                        eval(['header.',fieldname,'= fieldvalue;']);
                    end
                end
            end
        end
       
    end

end








