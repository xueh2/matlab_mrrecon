
function timeUsed = UpdateGadgetronOnRemoteAll(clean, branch, gt_branch, gtprep_branch)
%% update gadgetron on remote server
% timeUsed = UpdateGadgetronOnRemoteAll(clean, gt_branch, gtprep_branch)
% update barbados, denmark, palau, samoa, andorra

if(nargin<1)
    clean = 0;
end

if(nargin<2)
    branch = 'master';
end

if(nargin<3)
    gt_branch = 'master';
end

if(nargin<4)
    gtprep_branch = 'master';
end

tic
% UpdateGadgetronOnRemote('barbados', clean, branch, gt_branch, gtprep_branch)
UpdateGadgetronOnRemote('denmark', clean, branch, gt_branch, gtprep_branch)
UpdateGadgetronOnRemote('bermuda', clean, branch, gt_branch, gtprep_branch)
UpdateGadgetronOnRemote('gibraltar', clean, branch, gt_branch, gtprep_branch)
UpdateGadgetronOnRemote('palau', clean, branch, gt_branch, gtprep_branch)
UpdateGadgetronOnRemote('hongkong', clean, branch, gt_branch, gtprep_branch)
UpdateGadgetronOnRemote('samoa', clean, branch, gt_branch, gtprep_branch)
% UpdateGadgetronOnRemote('andorra', clean, branch, gt_branch, gtprep_branch)
timeUsed = toc;