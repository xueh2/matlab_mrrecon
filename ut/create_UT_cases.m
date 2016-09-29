function UTCases = create_UT_cases(home, pathStr, VDorVB, configFile, method)
% create the UTCases for all data in the home folder
% VDorVB :'VD' or 'VB'
% configFile: the gadgetron configuration file
% method: 'grappa', 'spirit' or 'l1spirit'

[cases, num] = FindSubDirs(home);

UTCases = cell(num, 6);

for ii=1:num
    UTCases{ii, 1} = pathStr;
    UTCases{ii, 2} = cases{ii};
    UTCases{ii, 3} = VDorVB;
    UTCases{ii, 4} = configFile;
    
    if ( strcmp(method, 'grappa') == 1 )
        UTCases{ii, 5} = 'grappa_res';
        UTCases{ii, 6} = 'grappa_ref';
    elseif ( strcmp(method, 'spirit') == 1 )
        UTCases{ii, 5} = 'spirit_res';
        UTCases{ii, 6} = 'spirit_ref';
    elseif ( strcmp(method, 'l1spirit') == 1 )    
        UTCases{ii, 5} = 'l1spirit_res';
        UTCases{ii, 6} = 'l1spirit_ref';
    end
end
     
UTCases
