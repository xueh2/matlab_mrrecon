
% Compute the moco performance statistics
% frame_index, dice, FP, FN, BSE
% frame_index_KLT, dice_KLT, FP_KLT, FN_KLT, BSE_KLT
% frame_index_progressive, dice_progressive, FP_progressive, FN_progressive, BSE_progressive

[mean_dice, max_dice, min_dice, std_dice] = ComputeMoCoPerformance_Perfusion_KLT_evaluation_mean_all(dice);
[mean_dice_KLT, max_dice_KLT, min_dice_KLT, std_dice_KLT] = ComputeMoCoPerformance_Perfusion_KLT_evaluation_mean_all(dice_KLT);
[mean_dice_progressive, max_dice_progressive, min_dice_progressive, std_dice_progressive] = ComputeMoCoPerformance_Perfusion_KLT_evaluation_mean_all(dice_progressive); 

min_dice_epi = zeros(3, 3);
min_dice_epi(:, 1) = min_dice(:, 1);
min_dice_epi(:, 2) = min_dice_progressive(:, 1);
min_dice_epi(:, 3) = min_dice_KLT(:, 1);

mean_dice_epi = zeros(3, 3);
mean_dice_epi(:, 1) = mean_dice(:, 1);
mean_dice_epi(:, 2) = mean_dice_progressive(:, 1);
mean_dice_epi(:, 3) = mean_dice_KLT(:, 1);

Labels = {'Basal', 'Apex', 'Medial'};

figure;
bar(min_dice_epi, 'grouped');
title('Minimal Dice');
legend('Ori', 'Direct', 'KLT');
set(gca, 'XTickLabel', Labels);
ylim([0.75 1])
 
figure;
bar(mean_dice_epi, 'grouped');
title('Mean Dice');
legend('Ori', 'Direct', 'KLT');
set(gca, 'XTickLabel', Labels);
ylim([0.75 1])

for s=1:3
    for e=1:2
        figure;
        hold on
        p = dice{s, e};
        ind = find(p(:)>0);
        
    end
end

blandAltmanPlot(seg_low_res(:,2), rt_low_res(:, 2), [110 270], [-25 20], FontSize, MarkerSize);
title('EDV (ml) 192*128');
set(gcf, 'OuterPosition', [300 400 512 360]);


[h, p] = ttest(seg_high_res_SAX(:,1), rt_high_res_SAX(:, 1))