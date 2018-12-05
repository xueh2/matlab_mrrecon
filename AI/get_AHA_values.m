
function [basal, mid, apex] = get_AHA_values(data, bulls_eye, max_d, min_d)
% get AHA mean/max/min/median/std for each segment

if(nargin<3)
    max_d = 1e5;
end
if(nargin<4)
    min_d = -1;
end

num_sectors = 1;
a = struct('m', zeros(num_sectors,1), 'median', zeros(num_sectors,1), 'max', zeros(num_sectors,1), 'min', zeros(num_sectors,1), 'sd', zeros(num_sectors,1), 'pixels', cell(num_sectors,1), ...
    'm_endo', zeros(num_sectors,1), 'median_endo', zeros(num_sectors,1), 'max_endo', zeros(num_sectors,1), 'min_endo', zeros(num_sectors,1), 'sd_endo', zeros(num_sectors,1), 'pixels_endo', [], ...
    'm_epi', zeros(num_sectors,1), 'median_epi', zeros(num_sectors,1), 'max_epi', zeros(num_sectors,1), 'min_epi', zeros(num_sectors,1), 'sd_epi', zeros(num_sectors,1), 'pixels_epi', []);

basal = [a;a;a;a;a;a];

mid = basal;

num_sectors = 1;
a = struct('m', zeros(num_sectors,1), 'median', zeros(num_sectors,1), 'max', zeros(num_sectors,1), 'min', zeros(num_sectors,1), 'sd', zeros(num_sectors,1), 'pixels', cell(num_sectors,1), ...
    'm_endo', zeros(num_sectors,1), 'median_endo', zeros(num_sectors,1), 'max_endo', zeros(num_sectors,1), 'min_endo', zeros(num_sectors,1), 'sd_endo', zeros(num_sectors,1), 'pixels_endo', [], ...
    'm_epi', zeros(num_sectors,1), 'median_epi', zeros(num_sectors,1), 'max_epi', zeros(num_sectors,1), 'min_epi', zeros(num_sectors,1), 'sd_epi', zeros(num_sectors,1), 'pixels_epi', []);

apex = [a;a;a;a];

for ii=1:3
    
    sectors = bulls_eye(ii).sectors;
    
    d = data(:,:, bulls_eye(ii).slice_index);    
    num_sectors = numel(bulls_eye(ii).contours);
        
    % 16 segments
    for s=1:num_sectors 
               
        ind = find(bulls_eye(ii).sectors(:)==s);
        pts = d(ind(:));        
        pts = apply_max_min(pts, max_d, min_d);
        
        if(ii==1)
            basal(s).m = mean(pts);
            basal(s).median = median(pts);
            basal(s).max = max(pts);
            basal(s).min = min(pts);
            basal(s).sd = std(pts);
            basal(s).pixels = pts;
        elseif (ii==2)
            mid(s).m = mean(pts);
            mid(s).median = median(pts);
            mid(s).max = max(pts);
            mid(s).min = min(pts);
            mid(s).sd = std(pts);
            mid(s).pixels = pts;
        else
            apex(s).m = mean(pts);
            apex(s).median = median(pts);
            apex(s).max = max(pts);
            apex(s).min = min(pts);
            apex(s).sd = std(pts);
            apex(s).pixels = pts;
        end
    end
    
    % 32 segments
    for s=1:num_sectors 
        
        pts = get_endo_epi_pts(d, bulls_eye(ii).endo_contours{s}, max_d, min_d);
                
        if(ii==1)
            basal(s).m_endo = mean(pts);
            basal(s).median_endo = median(pts);
            basal(s).max_endo = max(pts);
            basal(s).min_endo = min(pts);
            basal(s).sd_endo = std(pts);
            basal(s).pixels_endo = pts;
        elseif (ii==2)
            mid(s).m_endo = mean(pts);
            mid(s).median_endo = median(pts);
            mid(s).max_endo = max(pts);
            mid(s).min_endo = min(pts);
            mid(s).sd_endo = std(pts);
            mid(s).pixels_endo = pts;
        else
            apex(s).m_endo = mean(pts);
            apex(s).median_endo = median(pts);
            apex(s).max_endo = max(pts);
            apex(s).min_endo = min(pts);
            apex(s).sd_endo = std(pts);
            apex(s).pixels_endo = pts;
        end
    end
    
    for s=1:num_sectors 
        
        pts = get_endo_epi_pts(d, bulls_eye(ii).epi_contours{s}, max_d, min_d);
                
        if(ii==1)
            basal(s).m_epi = mean(pts);
            basal(s).median_epi = median(pts);
            basal(s).max_epi = max(pts);
            basal(s).min_epi = min(pts);
            basal(s).sd_epi = std(pts);
            basal(s).pixels_epi = pts;
        elseif (ii==2)
            mid(s).m_epi = mean(pts);
            mid(s).median_epi = median(pts);
            mid(s).max_epi = max(pts);
            mid(s).min_epi = min(pts);
            mid(s).sd_epi = std(pts);
            mid(s).pixels_epi = pts;
        else
            apex(s).m_epi = mean(pts);
            apex(s).median_epi = median(pts);
            apex(s).max_epi = max(pts);
            apex(s).min_epi = min(pts);
            apex(s).sd_epi = std(pts);
            apex(s).pixels_epi = pts;
        end
    end
end

end

function pts = apply_max_min(pts, max_d, min_d)
    ind2 = find(pts>=min_d & pts<=max_d);
    pts = pts(ind2);
end

function pts = get_endo_epi_pts(d, endo, max_d, min_d)
    BW=roipoly(d, endo(:,1), endo(:,2));
    ind=find(BW >0);        
    pts = d(ind(:));        
    pts = apply_max_min(pts, max_d, min_d);
end