
function [key, user] = sshKeyLookup(host)
%% find the ssh key for host
% [key, user] = sshKeyLookup(host)

key_dir = getenv('GADGETRON_KEY_FOLDER');

key = [];
user = [];

if ( strcmp(host, 'barbados') == 1)
    key = fullfile(key_dir, 'xueh2_barbados_private');
    user = 'xueh2';
elseif ( strcmp(host, 'denmark') == 1 )
    key = fullfile(key_dir, 'gtuser_denmark_private');
    user = 'gtuser';
elseif ( strcmp(host, 'palau') == 1 )
    key = fullfile(key_dir, 'xueh2_palau_private');
    user = 'xueh2';
elseif ( strcmp(host, 'nauru') == 1 )
    key = fullfile(key_dir, 'xueh2_nauru_private');
    user = 'xueh2';    
elseif ( strcmp(host, 'samoa') == 1 )
    key = fullfile(key_dir, 'xueh2_samoa_private');
    user = 'xueh2';
elseif ( strcmp(host, 'andorra') == 1 )
    key = fullfile(key_dir, 'xueh2_andorra_private');
    user = 'xueh2';  
elseif ( strcmp(host, 'hongkong') == 1 )
    key = fullfile(key_dir, 'xueh2_hongkong_private');
    user = 'xueh2';   
elseif ( strcmp(host, 'bermuda') == 1 )
    key = fullfile(key_dir, 'xueh2_bermuda_private');
    user = 'xueh2';       
elseif ( strcmp(host, 'gibraltar') == 1 )
    key = fullfile(key_dir, 'xueh2_gibraltar_private');
    user = 'xueh2';       
elseif ( strcmp(host, 'grenada') == 1 )
    key = fullfile(key_dir, 'xueh2_grenada_private');
    user = 'xueh2';       
elseif ( strcmp(host, 'localhost') == 1 )
    key = 'none';
    user = 'xueh2';
end
    