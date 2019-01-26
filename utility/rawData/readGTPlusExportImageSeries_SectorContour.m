
function [data, sector_c] = readGTPlusExportImageSeries_SectorContour(folderName, file_name, plot_flag)
% [data, sector_c] = readGTPlusExportImageSeries_SectorContour(folderName, file_name)

if(nargin<3)
    plot_flag = 0;
end

data = [];
sector_c = [];

real_name = fullfile(folderName, file_name);
data = analyze75read([real_name '.hdr']);

RO = size(data, 1);
E1 = size(data, 2);
            
xmlContent = xml_load(fullfile(folderName, [file_name '.attrib']));

N = numel(xmlContent);
for n=1:N
    if ( ~isempty(strfind(xmlContent(n).meta(1).name, 'GT_ROI')) )
        
        len = length('GT_ROI_S');
        
        sector_no = str2double(xmlContent(n).meta(1).name(len+1:end));                        
        
        num_pt = numel(xmlContent(n).meta)-4;
        curr_epi_pt = [];
        for pt=1:2:num_pt
            curr_epi_pt = [curr_epi_pt; str2double(xmlContent(n).meta(4+pt).value) str2double(xmlContent(n).meta(4+pt+1).value)];
        end
        
        sector_c = [sector_c; {sector_no, curr_epi_pt}];
    end
end

if(plot_flag)
    figure; 
    imagescn(permute(data, [2 1]));
    hold on
    for p=1:size(sector_c, 1)
        sn = sector_c{p, 1};
        sc = sector_c{p, 2};
        
        if(sn==1 | sn == 7 | sn == 13)
            plot(sc(:,2), sc(:,1), 'y-', 'LineWidth', 4);
        elseif(sn==2 | sn == 8 | sn == 14)
            plot(sc(:,2), sc(:,1), 'g-', 'LineWidth', 4);
        else
            plot(sc(:,2), sc(:,1), 'r-', 'LineWidth', 2);
        end
    end
    hold off
end