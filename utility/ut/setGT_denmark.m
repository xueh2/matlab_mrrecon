function setGT_denmark(port)
% setGT(server, port)

server = 'denmark'

if(nargin<1)
    port = 9002;
end

setenv('GT_HOST', server); setenv('GT_PORT', num2str(port));
getenv('GT_HOST')
getenv('GT_PORT')