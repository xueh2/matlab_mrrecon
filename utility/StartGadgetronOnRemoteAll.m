
function timeUsed = StartGadgetronOnRemoteAll()
%% start gadgetron on remote server
% timeUsed = StartGadgetronOnRemoteAll()

tic;

StartGadgetronOnRemote('barbados', 9008);
StartGadgetronOnRemote('denmark', 9008);
StartGadgetronOnRemote('bermuda', 9008);
StartGadgetronOnRemote('gibraltar', 9008);
StartGadgetronOnRemote('hongkong', 9008);
StartGadgetronOnRemote('palau', 9008);
StartGadgetronOnRemote('samoa', 9008);
StartGadgetronOnRemote('andorra', 9008);

timeUsed = toc;