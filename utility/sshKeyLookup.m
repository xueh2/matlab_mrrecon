
function [key, user] = sshKeyLookup(host)
%% find the ssh key for host
% [key, user] = sshKeyLookup(host)

key_dir = getenv('GADGETRON_KEY_FOLDER');

key = [];
user = [];

if isunix; user = getenv('USER');else; user = getenv('username');end

if ( strcmp(host, 'barbados') == 1 | strcmp(host, '137.187.135.191')==1)
    key = fullfile(key_dir, 'xueh2_barbados_private');
    user = 'xueh2';

elseif ( strcmp(host, 'denmark') == 1 )
    key = fullfile(key_dir, 'xueh2_private');

    
elseif ( strcmp(host, 'bahrain') == 1 | strcmp(host, '137.187.135.36')==1 )
    key = fullfile(key_dir, 'xueh2_private');

elseif ( strcmp(host, 'palau') == 1 )
    key = fullfile(key_dir, 'xueh2_palau_private');
 
elseif ( strcmp(host, 'nauru') == 1 )
    key = fullfile(key_dir, 'xueh2_nauru_private');
   
elseif ( strcmp(host, 'samoa') == 1 )
    key = fullfile(key_dir, 'xueh2_samoa_private');
 
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
       
elseif ( strcmp(host, 'grenada') == 1 )
    key = fullfile(key_dir, 'xueh2_grenada_private');
    user = 'xueh2';
elseif ( strcmp(host, 'localhost') == 1 )
    key = 'none';

elseif ( ~isempty(strfind(host, 'xueh2worker2.eastus.cloudapp.azure.com')) )    
    key = fullfile(key_dir, 'gtuser_gtCUDA_private');
    user = 'xueh2';
    
elseif ( ~isempty(strfind(host, 'eastus.cloudapp.azure.com')) )
    key = fullfile(key_dir, 'gtuser_gtCUDA_private');
    user = 'xueh2';
elseif ( ~isempty(strfind(host, 'localhost')) )
    key = '/home/xueh2/.ssh/id_rsa';
    user = 'xueh2';
elseif(strcmp(host, '137.187.135.157'))
    if(isunix())
        [~, result] = system('ifconfig');
    else
        [~, result] = system('ipconfig');
    end
    
    if(~isempty(strfind(result, '137.187.135.157')))
        key = '/home/boss/.ssh/id_rsa';
        user = 'boss';
    else
        key = '/home/xueh2/.ssh/id_rsa';
        user = 'xueh2';
    end
    
    if ~isempty(strfind(getenv('USER'),'kellmanp'))
        key = '/home/kellmanp/.ssh/beast_kellman';
        user = 'kellmanp';
    end
    
    
elseif ( ~isempty(strfind(host, '137.187.134.169')) | ~isempty(strfind(host, 'gt1')) | ~isempty(strfind(host, 'gt2') ) )
    key = '/home/xueh2/.ssh/id_rsa';
    user = 'xueh2';    
end

ind = find(key=='\');
if(~isempty(ind))
    key(ind) = '/';
end
