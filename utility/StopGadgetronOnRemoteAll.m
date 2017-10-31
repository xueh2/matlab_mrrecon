
function timeUsed = StopGadgetronOnRemoteAll()
%% stop gadgetron on remote server
% timeUsed = StopGadgetronOnRemoteAll()

tic;

% StopGadgetronOnRemote('barbados', 9008);
StopGadgetronOnRemote('denmark', 9008);
% StopGadgetronOnRemote('palau', 9008);
StopGadgetronOnRemote('samoa', 9017);
% StopGadgetronOnRemote('andorra', 9008);
StopGadgetronOnRemote('bermuda', 9008);
StopGadgetronOnRemote('hongkong', 9008);
StopGadgetronOnRemote('gibraltar', 9008);
StopGadgetronOnRemote('grenada', 9008);

timeUsed = toc;