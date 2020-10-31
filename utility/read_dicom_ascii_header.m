function header=read_dicom_ascii_header(filename, displaytextflag);
% function header=read_dicom_ascii_header(filename);
%
% function to read siemens ascii protocol header from dicom image file
% and output the information in a matlab structure variable "header"

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

% ascii file begins with:   ### ASCCONV BEGIN ###
% and ends with:            ### ASCCONV END ###

startflag=0;
% stopflag=0;
header=[];
start=[]; stop=[];
if nargin <2; displaytextflag = 0; end

fid=fopen(filename);
cnt=1;
while 1
    tline = fgetl(fid);
        tmp=char(tline);
%         start=findstr('### ASCCONV BEGIN ###',tmp);
%         stop=findstr('### ASCCONV END ###',tmp);
        start=strfind(tmp, '### ASCCONV BEGIN ###');
%         if ~isempty(start)
            stop=strfind(tmp, '### ASCCONV END ###');
%             if ~isempty(stop)
%                 disp(num2str(cnt));
%                 return
%             end
            if startflag==0
                stop=[];
            end
            
%         end
        % detect dicoms in step mode (overlays)
%         step=findstr('ENDSTEP;',tmp);
        step=strfind(tmp,'ENDSTEP;');
        if ~isempty(step);
            disp('found ENDSTEP;')
            fclose(fid); 
            return;
        end
        
    %     startflag=or(startflag,~isempty(start));
        if and(startflag,and(~isempty(stop),cnt>20));
%         if and(~isempty(stop),cnt>50);
%             ~isempty(start)
            fclose(fid); 
            return;
        end
        if startflag
        %     if ~ischar(tline); fclose(fid); return; end
            index=findstr(tline,'=');
            if displaytextflag==1
                disp(tline); end
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
                if ~isempty(index2)
                    fieldname = strrep(fieldname,fieldname(index2+1:index3-1),...
                        num2str(str2num(fieldname(index2+1:index3-1))+1));
                end
                try
                    eval(['header.',fieldname,'= fieldvalue;']);
                catch
                    start=[];
                end
            end
        end
        startflag=or(startflag,~isempty(start));
        cnt=cnt+1;
end
fclose(fid);
