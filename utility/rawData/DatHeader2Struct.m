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
