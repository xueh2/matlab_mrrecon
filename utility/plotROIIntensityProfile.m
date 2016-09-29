function plotROIIntensityProfile(Profile_info_table, lineWidth)

% N = size(Profile_info_table, 2);
% M = size(Profile_info_table, 1);
% 
% for m=1:M
%     figure;
%     hold on
% 
%     plot(Profile_info_table(m,1).Profile_Data_X, Profile_info_table(m,1).Profile_Data_Y, 'b-', 'LineWidth',lineWidth);
%     plot(Profile_info_table(m,2).Profile_Data_X, Profile_info_table(m,2).Profile_Data_Y, 'k:', 'LineWidth',lineWidth);
%     plot(Profile_info_table(m,3).Profile_Data_X, Profile_info_table(m,3).Profile_Data_Y, 'k-.', 'LineWidth',lineWidth);
%     plot(Profile_info_table(m,4).Profile_Data_X, Profile_info_table(m,4).Profile_Data_Y, 'k--', 'LineWidth',lineWidth);
%     plot(Profile_info_table(m,5).Profile_Data_X, Profile_info_table(m,5).Profile_Data_Y, 'r-', 'LineWidth',lineWidth);
%     
%     hold off
%     box on
%     % legend('Grappa', 'Grappa-MoCo', 'L1SPIRiT', 'NR-SPIRiT');
%     legend('GroundTruth', 'Grappa', 'Grappa-MoCo', 'L1SPIRiT', 'NR-SPIRiT');
%     title(['Profile ' num2str(m)]);
% 
% end

N = size(Profile_info_table, 2);
M = size(Profile_info_table, 1);

for n=1:N
    figure;
    hold on

    plot(Profile_info_table(1,n).Profile_Data_X, Profile_info_table(1,n).Profile_Data_Y, 'r', 'LineWidth',lineWidth);
    plot(Profile_info_table(2,n).Profile_Data_X, Profile_info_table(2,n).Profile_Data_Y, 'b', 'LineWidth',lineWidth);
    plot(Profile_info_table(3,n).Profile_Data_X, Profile_info_table(3,n).Profile_Data_Y, 'g', 'LineWidth',lineWidth);
    legend('LV', 'Myo', 'RV');
    
    hold off
    box on
    set(gca, 'FontSize', 14);
end
