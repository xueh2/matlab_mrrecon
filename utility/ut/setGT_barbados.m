function setGT_barbados(port)
% setGT(server, port)

server = 'barbados'

if(nargin<1)
    port = 9006;
end

setenv('GT_HOST', server); setenv('GT_PORT', num2str(port));
getenv('GT_HOST')
getenv('GT_PORT')