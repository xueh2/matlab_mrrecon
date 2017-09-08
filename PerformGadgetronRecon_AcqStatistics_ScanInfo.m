
function week_record = PerformGadgetronRecon_AcqStatistics_ScanInfo(scan_info, sites, fig_dir)
% week_record = PerformGadgetronRecon_AcqStatistics_ScanInfo(scan_info)

if(nargin<2)
    sites = [];
end

if(nargin < 3)
    fig_dir = '\\hl-share\RawMRI\Lab-Kellman\Share\data_record';
end

num = numel(scan_info);

patients = zeros(num, 1);
scans = zeros(num, 7); % LGE, DBLGE, Perf, Binning, RTCine, T2W, T2S

minw = zeros(num, 1);
maxw = zeros(num, 1);

first_date = 99999999;
last_date = 0;

for ii=1:num
    
    if(~isempty(sites))
        disp(['Scan info for ' sites{ii}]);
    end
    res = PerformGadgetronRecon_parse_scan_info(scan_info{ii});
    
    patients(ii) = res.num_pt;
    scans(ii, 1) = res.LGE;
    scans(ii, 2) = res.DBLGE;
    scans(ii, 3) = res.Perf;
    scans(ii, 4) = res.Binning;
    scans(ii, 5) = res.RTCine;
    scans(ii, 6) = res.T2W;
    scans(ii, 7) = res.T2S;     
    
    minw(ii) = res.week_num(1);
    maxw(ii) = res.week_num(2);
    
    if(res.study_dates(1)<first_date)
        first_date = res.study_dates(1);
    end
    
    if(res.study_dates(2)>last_date)
        last_date = res.study_dates(2);
    end
end

tScan = sum(scans, 1);
tScan = tScan(:);

disp(['Total number of patients : ' num2str(sum(patients))]);
disp(['Total number of LGE      : ' num2str(tScan(1))]);
disp(['Total number of DBLGE    : ' num2str(tScan(2))]);
disp(['Total number of Perf     : ' num2str(tScan(3))]);
disp(['Total number of Binning  : ' num2str(tScan(4))]);
disp(['Total number of RTCine   : ' num2str(tScan(5))]);
disp(['Total number of T2W      : ' num2str(tScan(6))]);
disp(['Total number of T2S      : ' num2str(tScan(7))]);

%% plot patients

labels = cell(num, 1);
for ii=1:num
    labels{ii} = [sites{ii} ', N=' num2str(patients(ii))];
end

h = figure('Name',['Number of Patients scanned with Gadgetron up to ' date],'NumberTitle','off');
pie(patients, labels);
title(['Number of Patients scanned with Gadgetron is ' num2str(sum(patients)) ', from ' num2str(first_date) ' to ' num2str(last_date)], 'FontSize', 16);
set(h, 'Position', [200 200 1.5*1024 1.5*768]);
saveas(h, fullfile(fig_dir, 'Patient_Count'), 'fig');

%% plot Perf

labels = cell(num, 1);
for ii=1:num
    labels{ii} = [sites{ii} ', N=' num2str(scans(ii, 3))];
end

h = figure('Name',['Number of Perfusion scans using Gadgetron up to ' date],'NumberTitle','off');
pie(scans(:, 3), labels);
title(['Number of Perfusion scans with Gadgetron is ' num2str(tScan(3)) ', from ' num2str(first_date) ' to ' num2str(last_date)], 'FontSize', 16);
set(h, 'Position', [300 200 1.5*1024 1.5*768]);
saveas(h, fullfile(fig_dir, 'Perfusion_Count'), 'fig');


%% plot LGE

labels = cell(num, 1);
for ii=1:num
    labels{ii} = [sites{ii} ', N=' num2str(scans(ii, 1))];
end

h = figure('Name',['Number of Bright blood LGE scans using Gadgetron up to ' date],'NumberTitle','off');
pie(scans(:, 1), labels);
title(['Number of Bright blood LGE scans with Gadgetron is ' num2str(tScan(1)) ', from ' num2str(first_date) ' to ' num2str(last_date)], 'FontSize', 16);
set(h, 'Position', [400 200 1.5*1024 1.5*768]);
saveas(h, fullfile(fig_dir, 'LGE_Count'), 'fig');

%% plot DB LGE

labels = cell(num, 1);
for ii=1:num
    labels{ii} = [sites{ii} ', N=' num2str(scans(ii, 2))];
end

h = figure('Name',['Number of Dark blood LGE scans using Gadgetron up to ' date],'NumberTitle','off');
pie(scans(:, 2), labels);
title(['Number of Dark blood LGE scans with Gadgetron is ' num2str(tScan(2)) ', from ' num2str(first_date) ' to ' num2str(last_date)], 'FontSize', 16);
set(h, 'Position', [500 200 1.5*1024 1.5*768]);
saveas(h, fullfile(fig_dir, 'DBLGE_Count'), 'fig');

%% count every week data

week_record = cell(num, 1);
for ii=1:num
    week_record{ii} = get_week_usage(scan_info{ii});
end

[bar_x, bar_y] = make_week_plot(week_record, min(minw), max(maxw), 'Perf');

h = figure('Name',['Number of perfusion scans using Gadgetron up to ' date],'NumberTitle','off');
hold on
bar(1:size(bar_y, 1), bar_y,'stacked','LineWidth',1.2);
set(gca,'XTick',1:size(bar_y, 1));
set(gca,'XTickLabel',bar_x)
set(gca,'FontSize', 14)
xlim([0 size(bar_y, 1)+1])
hold off
xlabel('Scan week num')
ylabel('# of Perfusion scans')
title(['Year 2016 - 2017, N=' num2str(sum(bar_y(:)))])
box on
legend(sites, 'Location', 'NorthWest' );
set(h, 'Position', [50 50 2.5*1024 1.5*768]);

saveas(h, fullfile(fig_dir, 'Perfusion_sites'), 'fig');
saveas(h, fullfile(fig_dir, 'Perfusion_sites'), 'tiff');

%%

[bar_x, bar_y] = make_week_plot(week_record, min(minw), max(maxw), 'LGE');

h = figure('Name',['Number of LGE scans using Gadgetron up to ' date],'NumberTitle','off');
hold on
bar(1:size(bar_y, 1), bar_y,'stacked','LineWidth',1.2);
set(gca,'XTick',1:size(bar_y, 1));
set(gca,'XTickLabel',bar_x)
set(gca,'FontSize', 14)
xlim([0 size(bar_y, 1)+1])
hold off
xlabel('Scan week num')
ylabel('# of LGE scans')
title(['Year 2016 - 2017, N=' num2str(sum(bar_y(:)))])
box on
legend(sites, 'Location', 'NorthWest' );
set(h, 'Position', [50 50 2.5*1024 1.5*768]);

saveas(h, fullfile(fig_dir, 'LGE_sites'), 'fig');
saveas(h, fullfile(fig_dir, 'LGE_sites'), 'tiff');

end

function week_record = get_week_usage(scan_info)

    num = size(scan_info, 1);
    
    weeks = [];
    
    for ii=1:num
        
        found = 0;
        for jj=1:numel(weeks)
            if(scan_info.week_num(ii)==weeks(jj))
                found = 1;
                break;
            end
        end
        
        if(~found)
            weeks = [weeks scan_info.week_num(ii)];
        end
    end

    week_record = zeros(numel(weeks), 7+1);
    week_record(:,1) = weeks(:);
    
    for ii=1:num        
        ind = find(scan_info.week_num(ii)==weeks);
        week_record(ind, 2) = week_record(ind, 2) + scan_info.LGE(ii);
        week_record(ind, 3) = week_record(ind, 3) + scan_info.DBLGE(ii);
        week_record(ind, 4) = week_record(ind, 4) + scan_info.Perf(ii);
        week_record(ind, 5) = week_record(ind, 5) + scan_info.Binning(ii);
        week_record(ind, 6) = week_record(ind, 6) + scan_info.RTCine(ii);
        week_record(ind, 7) = week_record(ind, 7) + scan_info.T2W(ii);
        week_record(ind, 8) = week_record(ind, 8) + scan_info.T2S(ii);
    end
end

function [bar_x, bar_y] = make_week_plot(week_record, minW, maxW, scan_to_plot)
% compute bar x and y for weekly plot for required scan

mW = weeknum(datenum('20161231', 'yyyymmdd'));

bar_x = [];
for kk=minW:maxW
    w = kk;
    if(w>mW) w = w - mW; end
    bar_x = [bar_x; {num2str(w)}];
end

num = numel(week_record);

bar_y = zeros(maxW-minW+1, num);

for ii=1:num    
    
    wr = week_record{ii};    
    numCase = size(wr, 1);
    
    for jj=1:numCase
        x = wr(jj, 1)-minW+1;
        
        if(strcmp(scan_to_plot, 'LGE')==1)
            bar_y(x, ii) = bar_y(x, ii) + wr(jj, 2);
            
        elseif(strcmp(scan_to_plot, 'DBLGE')==1)
            bar_y(x, ii) = bar_y(x, ii) + wr(jj, 3);
            
        elseif(strcmp(scan_to_plot, 'Perf')==1)
            bar_y(x, ii) = bar_y(x, ii) + wr(jj, 4);
            
        elseif(strcmp(scan_to_plot, 'Binning')==1)
            bar_y(x, ii) = bar_y(x, ii) + wr(jj, 5);
            
        elseif(strcmp(scan_to_plot, 'RTCine')==1)
            bar_y(x, ii) = bar_y(x, ii) + wr(jj, 6);
            
        elseif(strcmp(scan_to_plot, 'T2W')==1)
            bar_y(x, ii) = bar_y(x, ii) + wr(jj, 7);
            
        elseif(strcmp(scan_to_plot, 'T2S')==1)
            bar_y(x, ii) = bar_y(x, ii) + wr(jj, 8);
        end
        
    end
end

end
