function get_one_batch_file_AI(aiDir, script_name)
% get_one_batch_file_AI(aiDir, script_name)

[names, num] = findFILE(aiDir, [script_name '*']);

if(isunix())
    ext = '.sh';
else
    ext = '.bat';
end

fid = fopen(fullfile(aiDir, [script_name '_all' ext]), 'a+');

GT_HOME = getenv('GADGETRON_HOME');
fprintf(fid, '%s\n', ['cd ' GT_HOME '/share/gadgetron/python']);

for n=1:num
    fprintf(fid, '%s\n', names{n});
end