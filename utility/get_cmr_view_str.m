function view_str = get_cmr_view_str(proto_str)

proto_str = lower(proto_str);

view_str = 'UNK';
if(~isempty(strfind(proto_str, 'molli')))
    view_str = 'NON';
    return;
end

if(~isempty(strfind(proto_str, 'cine')))
    view_str = 'TBD';
end

if(~isempty(strfind(proto_str, '4ch')) | ~isempty(strfind(proto_str, 'ch4')) | ~isempty(strfind(proto_str, '4 ch')) | ~isempty(strfind(proto_str, '4 c')) | ~isempty(strfind(proto_str, 'ipat 4')))
    view_str = 'CH4';
else
    if (~isempty(strfind(proto_str, 'vla')) | ~isempty(strfind(proto_str, '2ch')) | ~isempty(strfind(proto_str, 'ch2')) | ~isempty(strfind(proto_str, '2 ch')) | ~isempty(strfind(proto_str, '2 c')) | ~isempty(strfind(proto_str, 'ipat 2')))
        view_str = 'CH2';
    else
        if (~isempty(strfind(proto_str, '3ch')) | ~isempty(strfind(proto_str, 'ch3')) | ~isempty(strfind(proto_str, '3 ch')) | ~isempty(strfind(proto_str, '3 c')) | ~isempty(strfind(proto_str, 'ipat 3')) )
            view_str = 'CH3';
        else
            if (~isempty(strfind(proto_str, 'sa')) | ~isempty(strfind(proto_str, 'sax')))
                view_str = 'SAX';
            end
        end
    end
end