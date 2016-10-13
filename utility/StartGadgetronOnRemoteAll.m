
function timeUsed = StartGadgetronOnRemoteAll()
%% start gadgetron on remote server
% timeUsed = StartGadgetronOnRemoteAll()

tic;

StartGadgetronOnRemote('barbados', 9008);
StartGadgetronOnRemote('denmark', 9008);
StartGadgetronOnRemote('palau', 9008);
StartGadgetronOnRemote('samoa', 9016);
StartGadgetronOnRemote('andorra', 9008);

timeUsed = toc;