
function [sf, rf, sf_i, rf_i] = PerformGadgetronRecon_SavedIsmrmrd_PerfStat_Summary(PerfTable, stress_column, rest_column, stress_flow_column, rest_flow_column, ischemia_flow_column, stress_aif_column, rest_aif_column, report_column)

sf = [];
rf = [];
sf_i = [];
rf_i = [];

s_aif = [];
s_aif_duration = [];
s_hb = [];

r_aif = [];
r_aif_duration = [];
r_hb = [];

pt_age = [];
pt_gender = [];

num_column = size(PerfTable, 2);
num = size(PerfTable, 1)-1;
for n=1:num
       
    stressCase = PerfTable{n+1, stress_column};
    restCase = PerfTable{n+1, rest_column};
    
    cmr_report = PerfTable{n+1, report_column};
    
    disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, stress_column} ' - ' PerfTable{n+1, rest_column}]); 
    disp(['==================================================================']);  
       
    if(~isnan(cmr_report) & ~isempty(strfind(cmr_report, 'reversed')))
        disp(['rest stress reversed !!!']); 
        stress_flow = PerfTable{n+1, rest_flow_column};
        rest_flow = PerfTable{n+1, stress_flow_column};
        ischemia_flow = PerfTable{n+1, ischemia_flow_column};

        sf = [sf; stress_flow];
        rf = [rf; rest_flow];
        sf_i = [sf_i; ischemia_flow];

        s_aif = [s_aif; PerfTable{n+1, rest_aif_column}];
        s_aif_duration = [s_aif_duration; PerfTable{n+1, 10}];
        s_hb = [s_hb; PerfTable{n+1, 11}];

        r_aif = [r_aif; PerfTable{n+1, stress_aif_column}];
        r_aif_duration = [r_aif_duration; PerfTable{n+1, 7}];
        r_hb = [r_hb; PerfTable{n+1, 8}];
    else
        stress_flow = PerfTable{n+1, stress_flow_column};
        rest_flow = PerfTable{n+1, rest_flow_column};
        ischemia_flow = PerfTable{n+1, ischemia_flow_column};

        sf = [sf; stress_flow];
        rf = [rf; rest_flow];
        sf_i = [sf_i; ischemia_flow];

        s_aif = [s_aif; PerfTable{n+1, stress_aif_column}];
        s_aif_duration = [s_aif_duration; PerfTable{n+1, 7}];
        s_hb = [s_hb; PerfTable{n+1, 8}];

        r_aif = [r_aif; PerfTable{n+1, rest_aif_column}];
        r_aif_duration = [r_aif_duration; PerfTable{n+1, 10}];
        r_hb = [r_hb; PerfTable{n+1, 11}];
    end
    
    pt_age = [pt_age; PerfTable{n+1, 12}];
    pt_gender = [pt_gender; PerfTable{n+1, 13}];

end

disp('=======================================================================');
disp(['Stress flow - ' num2str(mean(sf(:))) '+/-' num2str(std(sf(:)))]);
disp(['Rest flow - ' num2str(mean(rf(:))) '+/-' num2str(std(rf(:)))]);

disp(['Stress aif - ' num2str(mean(s_aif(:))) '+/-' num2str(std(s_aif(:)))]);
disp(['Rest aif - ' num2str(mean(r_aif(:))) '+/-' num2str(std(r_aif(:)))]);

disp(['Stress heart beat - ' num2str(mean(s_hb(:))) '+/-' num2str(std(s_hb(:)))]);
disp(['Rest heart beat - ' num2str(mean(r_hb(:))) '+/-' num2str(std(r_hb(:)))]);

ind = find(sf_i>0);
if(~isempty(ind))
    disp(['Stress flow, ischemia - ' num2str(mean(sf_i(ind))) '+/-' num2str(std(sf_i(ind)))]);
end
