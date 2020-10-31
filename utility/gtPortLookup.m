
function GT_PORT = gtPortLookup(gt_host)
%% find the port for gt host
% [gtPortLookup(gt_host)

GT_PORT = '9008';

if(strcmp(gt_host, 'palau'))
    GT_PORT = '9008';
end

if(strcmp(gt_host, 'localhost'))
    GT_PORT = '9002';
end

if(strcmp(gt_host, 'denmark'))
    GT_PORT = '9008';
end

if(strcmp(gt_host, 'samoa'))
    GT_PORT = '9016';
end

if(strcmp(gt_host, 'barbados') | strcmp(gt_host, '137.187.134.104'))
    GT_PORT = '9008';
end

if(strcmp(gt_host, 'andorra'))
    GT_PORT = '9008';
end

if(strcmp(gt_host, 'hongkong'))
    GT_PORT = '9008';
end

if(strcmp(gt_host, 'bermuda'))
    GT_PORT = '9008';
end

if(strcmp(gt_host, 'gibraltar'))
    GT_PORT = '9008';
end

if( ~isempty(strfind(gt_host, 'eastus.cloudapp.azure.com')) )
    GT_PORT = '9008';
end

if(strcmp(gt_host, '137.187.135.157'))
    if(isunix())
        [~, result] = system('ifconfig');
    else
        [~, result] = system('ipconfig');
    end
    
    if(~isempty(strfind(result, '137.187.135.157')))
        GT_PORT = '9006';
    else
        GT_PORT = '9008';
    end
end

if(strcmp(gt_host, 'gt1') | strcmp(gt_host, '137.187.135.169'))
    GT_PORT = '9008';
end


