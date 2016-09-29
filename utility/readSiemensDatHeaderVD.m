function [headers, prot] = readSiemensDatHeaderVD(h5filename, scanno)
% read in the siemens meas dat file header for VD

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
fprintf(fid, '%s', headers.Config);
fclose(fid);

fid = fopen(fullfile(path, [name '_Dicom.txt']), 'w');
fprintf(fid, '%s', headers.Dicom);
fclose(fid);

fid = fopen(fullfile(path, [name '_MeasYaps.txt']), 'w');
fprintf(fid, '%s', headers.MeasYaps);
fclose(fid);

fid = fopen(fullfile(path, [name '_Phoenix.txt']), 'w');
fprintf(fid, '%s', headers.Phoenix);
fclose(fid);

fid = fopen(fullfile(path, [name '_Spice.txt']), 'w');
fprintf(fid, '%s', headers.Spice);
fclose(fid);

fid = fopen(fullfile(path, [name '_MrcTime.txt']), 'w');
fprintf(fid, '%s', headers.MrcTime);
fclose(fid);

% prot = MeasYaps2structCall(headers);
prot = headers.MeasYaps;

end

% --------------------------------------------------------------

% header description
% // |X X X X|X X X X|X X X X X....0|X X X X|X X X X X.....|X X X X X X....0|X X X X|X X X X X.....|XXXX..............|XXXX....
% // |4 bytes|4 bytes|   x Bytes    |4 Bytes|     x Bytes  |    x Bytes     |4 Bytes|  x Bytes     | padding          | data
% // |hdr len|buf nr.| name 1       |len 1  | prot 1       | name 2         |len 2  | prot 2       | 32byte aligned   |
% 
% Description:
% First 4 Bytes: Overall length of header (hdr len).
% This can be used to hop directly to the raw data.
% The next 4 bytes indicates the number of embedded data structures (buf nr.)
% Each data structure starts with a name (e.g. name 1):
% this is a NULL terminated string, then the next 4 bytes are
% indicating the length (e.g. len 1) of the data structure and
% the data structure (e.g. prot 1) itself.
function header = MeasYaps2structCall(headers)
    % function header=MeasYaps2struct(headers);
    %
    % function to covert MeasYaps from text to Matlab structure

    %     ***************************************
    %     *  Peter Kellman  (kellman@nih.gov)   *
    %     *  Laboratory for Cardiac Energetics  *
    %     *  NIH NHLBI                          *
    %     ***************************************

    MeasYaps=headers.MeasYaps;
    lf = find(double(MeasYaps)==10); % find positions of line feeds
    lf = [0,lf];

    for line=1:length(lf)-1
        tline = MeasYaps(lf(line)+1:lf(line+1)-1);
        if ~ischar(tline); return; end
        index=findstr(tline,'=');

        if ~isempty(index)
            fieldname=deblank(tline(1:findstr(tline,'=')-1));
            fieldvalue=deblank(tline(findstr(tline,'=')+1:end));
            if max(isletter(fieldvalue))==0
                fieldvalue=str2num(fieldvalue);
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
            eval(['header.',fieldname,'= fieldvalue;']);
        end
    end
end
% --------------------------------------------------------------
