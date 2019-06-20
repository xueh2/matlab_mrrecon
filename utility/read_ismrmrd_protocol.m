function [patientID, protocol , ismrmrd_header] = read_ismrmrd_protocol(filename);
% function [patientID, protocol, ismrmrd_header] = read_ismrmrd_protocol(filename);

% uses Matlab Central base64decode.m
try
    dset = ismrmrd.Dataset(filename, 'dataset');
catch
    patientID = [];
    protocol = [];
    return
end
try
    hdr = ismrmrd.xml.deserialize(dset.readxml);
catch
   hdr.userParameters = []; 
end
if isfield(hdr.userParameters,'userParameterBase64')
    str = hdr.userParameters.userParameterBase64.value;
else
    patientID = [];
    protocol = [];
    return    
end

    
if isfield(hdr.subjectInformation,'patientID')
    patientID = hdr.subjectInformation.patientID;
else
    patientID = [];
end


% if nargout > 2
%     if isfield(hdr.userParameters,'userParameterDouble')
%         patinfo.PatientAge = hdr.userParameters.userParameterDouble(12);
%     else
%         patinfo.PatientAge = [];
%     end
%     if isfield(hdr.userParameters,'userParameterDouble')
%         patinfo.PatientWeight = hdr.userParameters.userParameterDouble(10);
%     else
%         patinfo.PatientWeight = [];
%     end
%     if isfield(hdr.userParameters,'userParameterDouble')
%         patinfo.PatientHeight = hdr.userParameters.userParameterDouble(11);
%     else
%         patinfo.PatientHeight = [];
%     end
% end
    
if nargout > 1
    header = char(base64decode(str));
    protocol = MeasYaps2struct(header);
end
if nargout > 2
    ismrmrd_header = hdr;
end


function protocol = MeasYaps2struct(header);
% function header=MeasYaps2struct(headers);
%
% function to covert MeasYaps from text to Matlab structure

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  NIH NHLBI                          *
%     ***************************************
warning off

MeasYaps = header;
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
protocol = header;
