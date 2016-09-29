function [Ki, Flow, Ki_2RR, Flow_2RR, Ki_NoR2Star, Flow_NoR2Star, Ki_NL, Flow_NL, Ki_2RR_NL, Flow_2RR_NL] = GetKiFlowValues(NL_ResultRecord)
% [Ki, Flow, Ki_2RR, Flow_2RR, Ki_NoR2Star, Flow_NoR2Star, Ki_NL, Flow_NL, Ki_2RR_NL, Flow_2RR_NL] = GetKiFlowValues(NL_ResultRecord)

N = size(NL_ResultRecord, 1)

% ----------------------

Ki = [];
Flow = [];

Ki_2RR = [];
Flow_2RR = [];

Ki_NoR2Star = [];
Flow_NoR2Star = [];

Ki_NL = [];
Flow_NL = [];

Ki_2RR_NL = [];
Flow_2RR_NL = [];

% ----------------------
% NL_ResultRecord = [NL_ResultRecord; {subdir{ii}, 'rest', Ki, Flow, nl_Ki, nl_Flow, half_Ki, half_Flow, nl_half_Ki, nl_half_Flow, Ki_without_R2Star, Flow_without_R2Star}];

for n=1:N

    % Ki
    v1 = NL_ResultRecord{n, 3};
    v2 = NL_ResultRecord{n, 5};
    v3 = NL_ResultRecord{n, 7};
    v4 = NL_ResultRecord{n, 9};
    v5 = NL_ResultRecord{n, 11};

    Ki = [Ki; v1];       
    Ki_NL = [Ki_NL; v2];
    Ki_2RR = [Ki_2RR; v3];
    Ki_2RR_NL = [Ki_2RR_NL; v3];
    Ki_NoR2Star = [Ki_NoR2Star; v5];  

    % Flow
    v1 = NL_ResultRecord{n, 4};
    v2 = NL_ResultRecord{n, 6};
    v3 = NL_ResultRecord{n, 8};
    v4 = NL_ResultRecord{n, 10};
    v5 = NL_ResultRecord{n, 12};

    Flow = [Flow; v1(:)];
    Flow_NL = [Flow_NL; v2(:)];
    Flow_2RR = [Flow_2RR; v3(:)];
    Flow_2RR_NL = [Flow_2RR_NL; v4(:)];
    Flow_NoR2Star = [Flow_NoR2Star; v5(:)];   
end
    