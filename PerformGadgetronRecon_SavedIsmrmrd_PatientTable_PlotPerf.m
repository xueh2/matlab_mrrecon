
function tUsed = PerformGadgetronRecon_SavedIsmrmrd_PatientTable_PlotPerf(PerfTable, resDir, onlyReview, ind)
% PerformGadgetronRecon_SavedIsmrmrd_PatientTable_PlotPerf(PerfTable, resDir, onlyReview)
% PerformGadgetronRecon_SavedIsmrmrd_PatientTable_PlotPerf(PerfTable, 'D:\data\ut\NewData\PaperResults\KAROLINSKA_Area_Aif_recorded', onlyReview)
% PerformGadgetronRecon_SavedIsmrmrd_PatientTable_PlotPerf(PerfTable, ''D:\data\ut\NewData\PaperResults\BARTS_Area_Aif_recorded', onlyReview)
% setenv('OutputFormat', 'h5')

if(nargin<4)
    ind = -1;
end

num = size(PerfTable, 1)-1;

stress_col = 28;
rest_col = 29;
flow_windowing = [0 6];

if(ind>0)
    closeall
    n = ind-1;
    disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, stress_col} ' - ' PerfTable{n+1, rest_col}]);       
    PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir,  PerfTable{n+1, stress_col}, PerfTable{n+1, rest_col}, flow_windowing, onlyReview);
else    
    tUsed = [];
    % for n=1:num
    for n=341:num
        disp([num2str(n) ' out of ' num2str(num) ' - Processing : ' PerfTable{n+1, stress_col} ' - ' PerfTable{n+1, rest_col}]);       
        PerformGadgetronRecon_Plot_PerfusionCase_StressRest(resDir,  PerfTable{n+1, stress_col}, PerfTable{n+1, rest_col}, flow_windowing, onlyReview);
        % if(onlyReview) pause; end
        pause
        closeall;
    end
end
