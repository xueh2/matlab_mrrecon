function setGT(server, port)
% setGT(server, port)

if(nargin<1)
    server = 'localhost'
end

if(nargin<2)
    port = 9002;
end

setenv('GT_HOST', server); setenv('GT_PORT', num2str(port));
getenv('GT_HOST')
getenv('GT_PORT')