
function PerformGadgetronRecon_AcqStatistics_PerfusionCase_StressRest(perf_cases)
% PerformGadgetronRecon_Statistics_PerfusionCase_StressRest(perf_cases)

num = size(perf_cases, 1);

patientID = [];
rest_time = [];
stress_time = [];
scannerID = [];

% count number of scans per study date
scan_count = [];

ind = 1;
for ii=1:num
    
    restCase = perf_cases{ii, 3}
    stressCase = perf_cases{ii, 2}
    
    [configName, scannerID_rest, patientID_rest, studyID_rest, measurementID_rest, study_dates, study_year, study_month, study_day, study_time_rest] = parseSavedISMRMRD(restCase);
    [configName, scannerID_stress, patientID_stress, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time_stress] = parseSavedISMRMRD(stressCase);
    
    patientID = [patientID; {patientID_stress}];
    scannerID = [scannerID; {scannerID_stress}];

    exist_date = 0;
    for jj=1:size(scan_count, 1)
        if(strcmp(scan_count{jj, 1}, study_dates)==1)
           exist_date = 1;
           break;
        end
    end
    
    if(exist_date)
        if(strcmp(scannerID_stress, '66016')==1)
           scan_count{jj, 3} = scan_count{jj, 3} + 1;
        else
            scan_count{jj, 4} = scan_count{jj, 4} + 1;
        end
    else
        if(strcmp(scannerID_stress, '66016')==1)
            scan_count = [scan_count; {study_dates scannerID_stress 1 0}];
        else
            scan_count = [scan_count; {study_dates scannerID_stress 0 1}];
        end
    end
end
 
disp(['Total number of stress/rest cases: ' num2str(num)]);

week_record = [];

for jj=1:size(scan_count, 1)
    sd = scan_count{jj,1};
    t = datenum(num2str(sd), 'yyyymmdd')
    weekn = weeknum(t);
    
    weekFound = 0;
    for kk=1:size(week_record, 1)
        if(week_record(kk,1)==weekn)
            weekFound = 1;
            break;
        end
    end
    
    if(weekFound)
        week_record(kk, 2) = week_record(kk, 2) + scan_count{jj,3};
        week_record(kk, 3) = week_record(kk, 3) + scan_count{jj,4};
    else
        week_record = [week_record; weekn scan_count{jj,3} scan_count{jj,4}];
    end
end

bar_x = [];
for jj=1:size(scan_count, 1)
    sd = scan_count{jj,1};
    bar_x = [bar_x; {sd(5:end)}];
end

bar_y = [];
for jj=1:size(scan_count, 1)
    bar_y = [bar_y; scan_count{jj,3} scan_count{jj,4}];
end

figure; 
hold on
bar(1:size(bar_y, 1), bar_y,'stacked', 'FaceColor',[0 .1 .7],'EdgeColor',[0.1 .1 .1],'LineWidth',0.2);
set(gca,'XTick',1:size(bar_y, 1));
set(gca,'XTickLabel',bar_x)
set(gca,'FontSize',6.0)
ylim([0 15])
xlim([1 size(bar_y, 1)])
hold off
xlabel('Scan date')
ylabel('# of stress/rest studies')
box on

%% by week
bar_x = [];
for kk=1:size(week_record, 1)
    bar_x = [bar_x; {num2str(week_record(kk,1))}];
end

bar_y = [];
for kk=1:size(week_record, 1)
    bar_y = [bar_y; week_record(kk,2) week_record(kk,3)];
end

figure; 
hold on
bar(1:size(bar_y, 1), bar_y,'stacked', 'FaceColor',[0 .1 .7],'EdgeColor',[0.1 .1 .1],'LineWidth',1.2);
set(gca,'XTick',1:size(bar_y, 1));
set(gca,'XTickLabel',bar_x)
xlim([0 size(bar_y, 1)+1])
hold off
xlabel('Scan week num')
ylabel('# of stress/rest studies')
title('Year 2016')
box on
