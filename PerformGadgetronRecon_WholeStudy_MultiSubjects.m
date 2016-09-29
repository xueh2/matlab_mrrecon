
function [tUsed, ignored] = PerformGadgetronRecon_WholeStudy_MultiSubjects(home, gt_host, startRemoteGT, deleteh5)
% perform gadgetron reconstruction for the whole study
% [tUsed, ignored] = PerformGadgetronRecon_WholeStudy_MultiSubjects(home, gt_host, startRemoteGT, deleteh5, start_case, end_case, ignore_list, styleSheet)

[names, num] = FindSubDirs(home);

tUsed = [];
ignored = [];

for n=1:num
    [tUsed_a, ignored_a] = PerformGadgetronRecon_WholeStudy(fullfile(home, names{n}), gt_host, startRemoteGT, deleteh5)
    
%     tUsed = [tUsed; tUsed_a];
%     ignored = [ignored; ignored_a];
end