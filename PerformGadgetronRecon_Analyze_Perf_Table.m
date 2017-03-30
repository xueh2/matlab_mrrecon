
function PerformGadgetronRecon_Analyze_Perf_Table(perf_tables, names)
% PerformGadgetronRecon_Analyze_Perf_Table(perf_tables, names)

N = numel(perf_tables);

color_all = ['r' 'b' 'k' 'm' 'g' 'y' 'c'];
symbol_all = ['+' 'o' '^' '*' 's' 'd' 'v'];

s_aif_peaks = cell(N,1);
s_aif_durations = cell(N,1);

r_aif_peaks = cell(N,1);
r_aif_durations = cell(N,1);

s_aif_peaks_3T = cell(N,1);
s_aif_durations_3T = cell(N,1);

r_aif_peaks_3T = cell(N,1);
r_aif_durations_3T = cell(N,1);

for ii=1:N
    
    perf_table = perf_tables{ii};
    name = names{ii};
    
    disp(['--> Analyze perf table from ' name]);
    disp(['================================================']);
    
    ind_case = perf_table.ind_case_all;
    ind_valid = find(ind_case>0);
    disp(['Found ' num2str(numel(ind_valid)) ' valid cases from ' num2str(numel(ind_case))]);
    
    field = perf_table.field_strength_all;
    ind_3T = find(field>2);
    disp(['Found ' num2str(numel(ind_3T)) ' 3T cases from ' num2str(numel(ind_valid))]);
    
    ind_1p5T = find(field>0.1 & field<2);
    
    disp(['---> 1.5T <---']);
    [s_aif_peak, s_aif_duration, r_aif_peak, r_aif_duration] = print_out_info(perf_table, ind_1p5T);    
    disp(['=============================']);
    disp(['---> 3T <---']);
    [s_aif_peak_3T, s_aif_duration_3T, r_aif_peak_3T, r_aif_duration_3T] = print_out_info(perf_table, ind_3T);    
    disp(['================================================']);    
    
    s_aif_peaks{ii} = s_aif_peak;
    s_aif_durations{ii} = s_aif_duration;
    
    r_aif_peaks{ii} = r_aif_peak;
    r_aif_durations{ii} = r_aif_duration;
    
    s_aif_peaks_3T{ii} = s_aif_peak_3T;
    s_aif_durations_3T{ii} = s_aif_duration_3T;
    
    r_aif_peaks_3T{ii} = r_aif_peak_3T;
    r_aif_durations_3T{ii} = r_aif_duration_3T;
end

figure
subplot(2,2,1)
plot_aif(names, s_aif_peaks, s_aif_durations, color_all, symbol_all, 'Stress AIF');
subplot(2,2,2)
plot_aif(names, r_aif_peaks, r_aif_durations, color_all, symbol_all, 'Rest AIF');
subplot(2,2,3)
plot_aif(names, s_aif_peaks_3T, s_aif_durations_3T, color_all, symbol_all, 'Stress AIF, 3T');
subplot(2,2,4)
plot_aif(names, r_aif_peaks_3T, r_aif_durations_3T, color_all, symbol_all, 'Rest AIF, 3T');

end

function [s_aif_peak, s_aif_duration, r_aif_peak, r_aif_duration] = print_out_info(perf_table, ind_1p5T)
    if(isempty(ind_1p5T))
        s_aif_peak = [];
        s_aif_duration = [];
        r_aif_peak = [];
        r_aif_duration = [];
        return;
    end
    
    v = perf_table.stress_aif_peak(ind_1p5T); ind = find(v<15); s_aif_peak = v(ind); v = v(ind);
    disp(['stress, aif peak Gd : ' num2str(mean(v)) ' +/- ' num2str(std(v))])
    v = perf_table.stress_aif_peak_no_T2Star(ind_1p5T);  v = v(ind);
    disp(['stress, aif peak Gd without T2* : ' num2str(mean(v)) ' +/- ' num2str(std(v))])
    v = perf_table.stress_aif_duration(ind_1p5T); v = v(ind); s_aif_duration = v;
    disp(['stress, aif first pass duration in seconds : ' num2str(mean(v)) ' +/- ' num2str(std(v))])
    v = perf_table.stress_aif_T2S_peak(ind_1p5T);  v = v(ind);
    disp(['stress, aif T2* at peak Gd : ' num2str(mean(v)) ' +/- ' num2str(std(v))])
    v = perf_table.HeartRate_stress(ind_1p5T);  v = v(ind);
    disp(['stress heart rate : ' num2str(mean(v)) ' +/- ' num2str(std(v))])
    disp(['--------------------']);
    v = perf_table.rest_aif_peak(ind_1p5T); ind = find(v<15); r_aif_peak = v(ind); v = v(ind);
    disp(['rest, aif peak Gd : ' num2str(mean(v)) ' +/- ' num2str(std(v))])
    v = perf_table.rest_aif_peak_no_T2Star(ind_1p5T); v = v(ind);
    disp(['rest, aif peak Gd without T2* : ' num2str(mean(v)) ' +/- ' num2str(std(v))])
    v = perf_table.rest_aif_duration(ind_1p5T); v = v(ind); r_aif_duration = v;
    disp(['rest, aif first pass duration in seconds : ' num2str(mean(v)) ' +/- ' num2str(std(v))])
    v = perf_table.rest_aif_T2S_peak(ind_1p5T); v = v(ind);
    disp(['rest, aif T2* at peak Gd : ' num2str(mean(v)) ' +/- ' num2str(std(v))])
    v = perf_table.HeartRate_rest(ind_1p5T); v = v(ind);
    disp(['rest heart rate : ' num2str(mean(v)) ' +/- ' num2str(std(v))])
    disp(['--------------------']);
    
    num = numel(ind_1p5T);
    
    r_s_time_diff = zeros(num, 1);
    for n=1:num
        st = perf_table.stress_time(ind_1p5T(n));
        st = st{1};
        rt = perf_table.rest_time(ind_1p5T(n));
        rt = rt{1};
        
        r_s_time_diff(n) = etime([0 0 0 str2num(rt(1:2)) str2num(rt(3:4)) str2num(rt(5:6))], [0 0 0 str2num(st(1:2)) str2num(st(3:4)) str2num(st(5:6))]);
    end
    
    v = r_s_time_diff;
    disp(['rest - stress scan time difference in seconds : ' num2str(mean(v)) ' +/- ' num2str(std(v))])
end

function plot_aif(names, aif_peak, aif_duration, color, symbol, title_str)

    names_used = [];

    % figure
    hold on
    
    for ii=1:numel(aif_peak)
        if(isempty(aif_peak{ii}))
            continue;
        end
        
        names_used = [names_used; names(ii)];
        plot(aif_duration{ii}, aif_peak{ii}, [color(ii) symbol(ii)], 'MarkerSize', 8);
    end    

    legend(names_used);
    
    for ii=1:numel(aif_peak)
        if(isempty(aif_peak{ii}))
            continue;
        end
                
        dd = [aif_duration{ii} aif_peak{ii}];
        m = mean(dd, 1);
        C = cov(dd);
        hp = plot_gaussian_ellipsoid(m, C, 1, 300, gca);
        set(hp, 'color', color(ii), 'LineWidth', 4); 
    end   
    hold off    
    xlabel('AIF duration (s)')
    ylabel('AIF peak Gd (mmol/L)');
    title(title_str)
    box on
    grid on
end