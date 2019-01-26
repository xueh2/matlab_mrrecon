
function [endoC, epiC, rvC, rviC] = perf_convert_whole_im_landmarks_to_scaled_cropped(endo, epi, rv, rvi, roi, rs)
% [endoC, epiC, rvC, rviC] = perf_convert_whole_im_landmarks_to_scaled_cropped(endo, epi, rv, rvi, roi, rs)
% input endo/epi/rv contours are 1 based cell array
% rvi: [SLC 2] rv insertion point, 1 based
% roi: [ps_x pe_x; ps_y pe_y; aif_s aif_e]
% rs: scale ratio

    endoC = endo;
    for n=1:numel(endo)
        endoC{n} = convert_contour(endo{n}, roi, rs);
    end
    
    epiC = epi;
    for n=1:numel(epi)
        epiC{n} = convert_contour(epiC{n}, roi, rs);
    end
    
    rvC = rv;
    for n=1:numel(rv)
        rvC{n} = convert_contour(rvC{n}, roi, rs);
    end
    
    rviC = rvi;
    rviC(1,:) = convert_contour(rvi(1,:), roi, rs);
    rviC(2,:) = convert_contour(rvi(2,:), roi, rs);
    rviC(3,:) = convert_contour(rvi(3,:), roi, rs);    
end

function C2 = convert_contour(C, roi, rs)

    if(isempty(C))
        C2 = [];
        return;
    end

    ps_x = roi(1,1);
    pe_x = roi(1,2);
    ps_y = roi(2,1);
    pe_y = roi(2,2);

    % rsC = rs*(C-1) + 1;
    rsC = rs*(C-1) + 1;
    
    rsC(:,1) = rsC(:,1) - ps_y + 1;
    rsC(:,2) = rsC(:,2) - ps_x + 1;
        
    C2 = rsC;
end
