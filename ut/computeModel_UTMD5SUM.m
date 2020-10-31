
function computeModel_UTMD5SUM(model_folder, search_str)
%% compute md5 sum for model files
% computeModel_UTMD5SUM(model_folder)

GTHome = getenv('GADGETRON_HOME')
UTDir = getenv('GTPLUS_UT_DIR')
OutputFormat = getenv('OutputFormat');

[m_files, num] = findFILE(model_folder, search_str);

for k=1:num

    disp(['#=====================================================================================']);
   
    ind = find(m_files{k}=='\');
    fname = m_files{k};
    fname(ind) = '/';
    % dos(['fciv.exe -md5 ' fname], '-echo');
    [status, outputs] = system(['fciv.exe -md5 ' fname]);
    print_json_record(outputs, fname);      
end

end

function print_json_record(outputs, fname)
    ind = strfind(outputs, '//');
    a = outputs(ind(end)+3:end);
    ind2 = strfind(a, ' ');
    sha1 = a(1:ind2(1)-1);
    [path, fname_str, ext] = fileparts(fname);
    disp(['{']);
    disp(['    "file":"' fname_str ext '",']);
    disp(['    "md5":"' sha1 '"']);
    disp(['}']);
end
