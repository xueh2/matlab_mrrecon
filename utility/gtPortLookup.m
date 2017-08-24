
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
    GT_PORT = '9008';
end

if(strcmp(gt_host, 'barbados'))
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
