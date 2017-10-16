
function timeUsed = StartGadgetronNHLBICloud()
%% start gadgetron on remote server for cloud
% timeUsed = StartGadgetronNHLBICloud()

tic;

% StopGadgetronOnRemote('barbados', 9008);
StopGadgetronOnRemote('denmark', 9008);
StopGadgetronOnRemote('bermuda', 9008);
StopGadgetronOnRemote('gibraltar', 9008);
StopGadgetronOnRemote('hongkong', 9008);
StopGadgetronOnRemote('palau', 9008);
StopGadgetronOnRemote('samoa', 9008);
StopGadgetronOnRemote('grenada', 9008);
% StopGadgetronOnRemote('andorra', 9008);

timeUsed = toc;

tic;

% StartGadgetronOnRemote('barbados', 9008);
StartGadgetronOnRemote('denmark', 9008);
StartGadgetronOnRemote('bermuda', 9008);
StartGadgetronOnRemote('gibraltar', 9008);
StartGadgetronOnRemote('hongkong', 9008);
% StartGadgetronOnRemote('palau', 9008);
StartGadgetronOnRemote('samoa', 9017);
StartGadgetronOnRemote('grenada', 9008);
% StartGadgetronOnRemote('andorra', 9008);

timeUsed = toc;

dos('gadgetron_cloudbus_relay &');