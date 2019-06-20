
function is_remote_computer = IsRemoteComputer(host)
% is_remote_computer = IsRemoteComputer(host)

if(isunix())
    [~, result] = system('ifconfig');
else
    [~, result] = system('ipconfig');
end 

is_remote_computer = 1;
if(~isempty(strfind(result, host)) | ~isempty(strfind(host, 'localhost')))
    is_remote_computer = 0;
end