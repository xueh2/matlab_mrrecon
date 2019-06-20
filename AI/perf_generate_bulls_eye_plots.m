
function [h_f, h_PS, h_Vb, h_Visf, h_Tc, h_delay, res] = perf_generate_bulls_eye_plots(bulls_eye, trainingCaseDir)
% create bulls eye plot for flow, PS, Visf, Vb, Tc and delay_time
% [h_f, h_PS, h_Vb, h_Visf, h_Tc, h_delay, res] = perf_generate_bulls_eye_plots(bulls_eye, trainingCaseDir)

% load data

f = readNPY(fullfile(trainingCaseDir, 'fmap_resized_training_norm.npy'));
PS = readNPY(fullfile(trainingCaseDir, 'PS_map_resized_training_norm.npy'));
Vb = readNPY(fullfile(trainingCaseDir, 'vb_map_resized_training_norm.npy'));
Visf = readNPY(fullfile(trainingCaseDir, 'visf_map_resized_training_norm.npy'));
Tc = readNPY(fullfile(trainingCaseDir, 'Tc_map_resized_training_norm.npy'));
delay = readNPY(fullfile(trainingCaseDir, 'delay_map_resized_training_norm.npy'));

RO = size(f, 1);
E1 = size(f, 2);

sRO = size(bulls_eye(1).sectors, 1);
sE1 = size(bulls_eye(1).sectors, 2);

if(RO~=sRO | E1~=sE1)
    f = Matlab_gt_resize_2D_image(double(f), sRO, sE1, 5);
    PS = Matlab_gt_resize_2D_image(double(PS), sRO, sE1, 5);
    Vb = Matlab_gt_resize_2D_image(double(Vb), sRO, sE1, 5);
    Visf = Matlab_gt_resize_2D_image(double(Visf), sRO, sE1, 5);
    Tc = Matlab_gt_resize_2D_image(double(Tc), sRO, sE1, 5);
    delay = Matlab_gt_resize_2D_image(double(delay), sRO, sE1, 5);
end

res = struct('basal_f', [], 'mid_f', [], 'apex_f', [], ...
    'basal_Vb', [], 'mid_Vb', [], 'apex_Vb', [], ...
    'basal_Visf', [], 'mid_Visf', [], 'apex_Visf', [], ...
    'basal_Tc', [], 'mid_Tc', [], 'apex_Tc', [], ...
    'basal_delay', [], 'mid_delay', [], 'apex_delay', [], ...
    'basal_PS', [], 'mid_PS', [], 'apex_PS', []);

cmap = PerfColorMap(0);
[h_f, h_f_32, res.basal_f, res.mid_f, res.apex_f] = create_bulls_eye(f, bulls_eye, [0 8], 'Flow, ml/min/g', cmap, 6, 0.1);

cmap = PSColorMap(0);
[h_PS, h_PS_32, res.basal_PS, res.mid_PS, res.apex_PS] = create_bulls_eye(PS, bulls_eye, [0 4], 'PS, ml/min/g', cmap, 4.5, 0.1);

cmap = MBVColorMap(0);
[h_Vb, h_Vb_32, res.basal_Vb, res.mid_Vb, res.apex_Vb] = create_bulls_eye(Vb, bulls_eye, [0 20], 'Vb, 100ml/g', cmap, 40, 1);

cmap = ECVColorMap(0);
[h_Visf, h_Visf_32, res.basal_Visf, res.mid_Visf, res.apex_Visf] = create_bulls_eye(Visf, bulls_eye, [0 80], 'Visf, 100ml/g', cmap, 80, 2);

cmap = PerfColorMap(0);
[h_Tc, h_Tc_32, res.basal_Tc, res.mid_Tc, res.apex_Tc] = create_bulls_eye(Tc, bulls_eye, [0 20], 'Tc, seconds', cmap, 25, 0.25);

cmap = PerfColorMap(0);
[h_delay, h_delay_32, res.basal_delay, res.mid_delay, res.apex_delay] = create_bulls_eye(delay, bulls_eye, [0 8], 'Delay, seconds', cmap, 40, 0.01);

end

function [h_16, h_32, basal, mid, apex] = create_bulls_eye(data, bulls_eye, windows, comments, cmap, max_d, min_d)

    [basal, mid, apex] = get_AHA_values(data, bulls_eye, max_d, min_d);
   
    h_16 = figure('Name', ['AHA 16, ' comments], 'NumberTitle','off', 'Renderer', 'painters', 'Position', [700 -800 1600 1600]);
    subplot(2,2,1)
    c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);
    set(c,'Color','k','LineWidth',2)

    fillBullseye(double(cat(2, apex.m)),0.5,1,-315,45); % inner circle (apical)
    fillBullseye(double(cat(2, mid.m)),1,1.5,-300,60); % mid circle (mid)
    fillBullseye(double(cat(2, basal.m)),1.5,2,-300,60); % outer circle (basal)

    uistack(c,'top');
    set(gca,'clim',windows);
    colormap(cmap);
    colorbar('South')
    axis off
    axis equal  
    title(comments, 'FontSize', 18)
    
    d = {'basal', basal(1).m, basal(2).m, basal(3).m, basal(4).m, basal(5).m, basal(6).m; ...
        'mid', mid(1).m, mid(2).m, mid(3).m, mid(4).m, mid(5).m, mid(6).m; ...
        'apex', apex(1).m, apex(2).m, apex(3).m, apex(4).m, 0, 0};
        
    uit = uitable(h_16, 'Data', d,'Position',[100 600 750 120], 'ColumnName', {'Slice, Mean', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6'}, 'FontSize', 16);
    
    d = {'basal', basal(1).sd, basal(2).sd, basal(3).sd, basal(4).sd, basal(5).sd, basal(6).sd; ...
        'mid', mid(1).sd, mid(2).sd, mid(3).sd, mid(4).sd, mid(5).sd, mid(6).sd; ...
        'apex', apex(1).sd, apex(2).sd, apex(3).sd, apex(4).sd, 0, 0};
        
    uit = uitable(h_16, 'Data', d,'Position',[100 450 750 120], 'ColumnName', {'Slice, SD', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6'}, 'FontSize', 16);
    
    d = {'basal', basal(1).max, basal(2).max, basal(3).max, basal(4).max, basal(5).max, basal(6).max; ...
        'mid', mid(1).max, mid(2).max, mid(3).max, mid(4).max, mid(5).max, mid(6).max; ...
        'apex', apex(1).max, apex(2).max, apex(3).max, apex(4).max, 0, 0};
        
    uit = uitable(h_16, 'Data', d,'Position',[100 300 750 120], 'ColumnName', {'Slice, Max', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6'}, 'FontSize', 16);
    
    d = {'basal', basal(1).min, basal(2).min, basal(3).min, basal(4).min, basal(5).min, basal(6).min; ...
        'mid', mid(1).min, mid(2).min, mid(3).min, mid(4).min, mid(5).min, mid(6).min; ...
        'apex', apex(1).min, apex(2).min, apex(3).min, apex(4).min, 0, 0};
        
    uit = uitable(h_16, 'Data', d,'Position',[100 150 750 120], 'ColumnName', {'Slice, Min', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6'}, 'FontSize', 16);
    
    subplot(2,2,2)
    
    h_32 = h_16;
%     h_32 = figure('Name', ['AHA 32, ' comments], 'NumberTitle','off');
    c1 = createBullseye([0 0.5 1 0; 0.5 0.75 4 45; 0.75 1 4 45; 1 1.25 6 0; 1.25 1.5 6 0; 1.5 1.75 6 0; 1.75 2 6 0],':');
    c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);
    
    set(c1,'Color','k','LineWidth',2)
    set(c,'Color','k','LineWidth',2)

    fillBullseye(double(cat(2, apex.m_endo)),0.5,0.75,-315,45); % inner circle (apical)
    fillBullseye(double(cat(2, mid.m_endo)),1.0,1.25,-300,60); % mid circle (mid)
    fillBullseye(double(cat(2, basal.m_endo)),1.5,1.75,-300,60); % outer circle (basal)
    
    fillBullseye(double(cat(2, apex.m_epi)),0.75,1,-315,45); % inner circle (apical)
    fillBullseye(double(cat(2, mid.m_epi)),1.25,1.5,-300,60); % mid circle (mid)
    fillBullseye(double(cat(2, basal.m_epi)),1.75,2,-300,60); % outer circle (basal)

%     c1 = createBullseye([0 0.5 1 0; 0.5 0.75 4 45; 0.75 1 4 45; 1 1.25 6 0; 1.25 1.5 6 0; 1.5 1.75 6 0; 1.75 2 6 0],':');
%     set(c1,'Color','k','LineWidth',2)
%     c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);
%     
%     set(c,'Color','k','LineWidth',2)
    
    uistack(c1,'top');   
    uistack(c,'top');
    set(gca,'clim',windows);
    colormap(cmap);
    colorbar('South')
    axis off
    axis equal
    title(comments, 'FontSize', 18)
    
    d = {'basal', basal(1).m_endo, basal(2).m_endo, basal(3).m_endo, basal(4).m_endo, basal(5).m_endo, basal(6).m_endo; ...
        'mid', mid(1).m_endo, mid(2).m_endo, mid(3).m_endo, mid(4).m_endo, mid(5).m_endo, mid(6).m_endo; ...
        'apex', apex(1).m_endo, apex(2).m_endo, apex(3).m_endo, apex(4).m_endo, 0, 0};
        
    uit = uitable(h_16, 'Data', d,'Position',[900 600 750 120], 'ColumnName', {'Slice, Mean, ENDO', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6'}, 'FontSize', 16);
    
    d = {'basal', basal(1).sd_endo, basal(2).sd_endo, basal(3).sd_endo, basal(4).sd_endo, basal(5).sd_endo, basal(6).sd_endo; ...
        'mid', mid(1).sd_endo, mid(2).sd_endo, mid(3).sd_endo, mid(4).sd_endo, mid(5).sd_endo, mid(6).sd_endo; ...
        'apex', apex(1).sd_endo, apex(2).sd_endo, apex(3).sd_endo, apex(4).sd_endo, 0, 0};
        
    uit = uitable(h_16, 'Data', d,'Position',[900 450 750 120], 'ColumnName', {'Slice, SD, ENDO', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6'}, 'FontSize', 16);
    
    d = {'basal', basal(1).m_epi, basal(2).m_epi, basal(3).m_epi, basal(4).m_epi, basal(5).m_epi, basal(6).m_epi; ...
        'mid', mid(1).m_epi, mid(2).m_epi, mid(3).m_epi, mid(4).m_epi, mid(5).m_epi, mid(6).m_epi; ...
        'apex', apex(1).m_epi, apex(2).m_epi, apex(3).m_epi, apex(4).m_epi, 0, 0};
        
    uit = uitable(h_16, 'Data', d,'Position',[900 300 750 120], 'ColumnName', {'Slice, Mean, EPI', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6'}, 'FontSize', 16);
    
    d = {'basal', basal(1).sd_epi, basal(2).sd_epi, basal(3).sd_epi, basal(4).sd_epi, basal(5).sd_epi, basal(6).sd_epi; ...
        'mid', mid(1).sd_epi, mid(2).sd_epi, mid(3).sd_epi, mid(4).sd_epi, mid(5).sd_epi, mid(6).sd_epi; ...
        'apex', apex(1).sd_epi, apex(2).sd_epi, apex(3).sd_epi, apex(4).sd_epi, 0, 0};
        
    uit = uitable(h_16, 'Data', d,'Position',[900 150 750 120], 'ColumnName', {'Slice, SD, EPI', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6'}, 'FontSize', 16);
end