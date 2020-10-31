

cd(home)
save parameters

% program setting

% ============================================================ %
% Global_Kmeans_MRF_Local_Kmeans_MRF_5classes
% clear
% load parameters

% Global_flag = 1; 
% Local_flag = 1; 
% 
% Atlas_flag_Global = 0;
% Atlas_flag_Local = 0;
% 
% PVs_flag = 1;
% Save_flag = 1;

% for i = 1:2
%     if ( i==1) 
%         Kmeans_flag_Global = 1;
%     end
%     if ( i==2) 
%         Kmeans_flag_Global = 0;
%     end
%     
%     for j = 1:2
%         if ( j==1 )
%             Kmeans_flag_Local = 1;
%         end
%         if ( j==2 )
%             Kmeans_flag_Local = 0;
%         end
        
% for k = 1:2
%     if (k==1)
%         MRF_flag_Global = 0;
%     end
%     if (k==2)
%         MRF_flag_Global = 1;
%     end
% 
%     for m = 1:2
%         if (m==1)
%             GMM_PVs_flag = 1;
%         end
%         if (m==2)
%             GMM_PVs_flag = 0;
%         end

        for n = 1:2
            if (n==1)
                Four_classes_flag = 1;
                Five_classes_flag = 0;
            end
            if (n==2)
                Four_classes_flag = 0;
                Five_classes_flag = 1;
            end
            
            Global_flag = 1; 
            Local_flag = 1;

            Atlas_flag_Global = 0;
            Atlas_flag_Local = 0;

            Kmeans_flag_Global = 1;
            Kmeans_flag_Local = 1;
            
            MRF_flag_Global = 0;
            MRF_flag_Local = 0;
            
            GMM_PVs_flag = 1;
            GMM_PVs_flag_Local = 1;
            
            Save_flag = 1;

            save settings
            
            clear
            load parameters
            load settings
            Neonatal_Brain_Segmentation
        end
%     end
% end

clear
load parameters