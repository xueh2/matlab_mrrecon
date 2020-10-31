
% ==================================================================== %
% Neonatal Brain Segmentation Script
% for the test of each step of current algorithm
% ==================================================================== %

% ==================================================================== %

% ==================================================================== %

% ==================================================================== %
% parse the parameters
prefix = CreatePrefix(home, Global_flag, Local_flag, Atlas_flag_Global, Atlas_flag_Local,...
    Kmeans_flag_Global, Kmeans_flag_Local, MRF_flag_Global, MRF_flag_Local, Four_classes_flag, ...
    Five_classes_flag, partsNumber, GMM_PVs_flag, GMM_PVs_flag_Local);
prefix
if ( Four_classes_flag )
    classnumber = 4;
end

if ( Five_classes_flag )
    classnumber = 5;
end

% ==================================================================== %
% load images

[brainmask, header] = LoadAnalyze(brainMaskfile,'Grey');
[imagedata, header] = LoadAnalyze(imagefile,'Grey');

% read kmeans intialization centres
if ( classnumber == 5 )
    Kmeansfile = [KmeansDir 'Kmeans_InitialCentres_5classes.txt'];
    KmeansResultFile = [KmeansDir '5classes_kmeansResult.hdr'];
end

if ( classnumber == 4 )
    Kmeansfile = [KmeansDir 'Kmeans_InitialCentres_4classes.txt'];
    KmeansResultFile = [KmeansDir '4classes_kmeansResult.hdr'];
end

% Kmeans_InitialCentres = ReadKmeans_InitialCentres(Kmeansfile, classnumber);
% Kmeans_InitialCentres
% ==================================================================== %
if ( Global_flag & ~Local_flag)
    % ============================================================== %
    if ( Kmeans_flag_Global )
             
        % detect PVs                          
        if ( classnumber == 5 )
            % run segmentation scripts
            if ( MRF_flag_Global )
                if ( GMM_PVs_flag )
                    Global_Kmeans_EM_MRF_5classes_GmmPVs;
                else
                    Global_Kmeans_EM_MRF_5classes;
                end
            else
                if ( GMM_PVs_flag )
                    Global_Kmeans_EM_5classes_PVs;
                else
                    Global_Kmeans_EM_5classes;
                end
            end
        end
        
        if ( classnumber == 4 )
            % run segmentation scripts
            if ( MRF_flag_Global )
                if ( GMM_PVs_flag )
                    Global_Kmeans_EM_MRF_4classes_GmmPVs;
                else
                    Global_Kmeans_EM_MRF_4classes;
                end
%                 Global_Kmeans_EM_MRF_4classes;
            else
                if ( GMM_PVs_flag )
                    Global_Kmeans_EM_4classes_PVs;
                else
                    Global_Kmeans_EM_4classes;
                end
            end
        end
    end
    % ============================================================== %
    if ( Atlas_flag_Global )
        disp('Atlas intialization has not been implemented yet ...');
    end
    % ============================================================== %
end

if ( ~Global_flag & Local_flag )
    % ============================================================== %
    if ( Kmeans_flag_Local )
             
        % detect PVs                          
        if ( classnumber == 5 )
            % run segmentation scripts
            if ( MRF_flag_Local )
                onlyLocal_Kmeans_EM_MRF_5classes;
            else
                onlyLocal_Kmeans_EM_5classes;
            end
        end
        
        if ( classnumber == 4 )
            % run segmentation scripts
            if ( MRF_flag_Local )
                onlyLocal_Kmeans_EM_MRF_4classes;
            else
                onlyLocal_Kmeans_EM_4classes;
            end
        end
    end
    % ============================================================== %
    if ( Atlas_flag_Local )
        disp('Atlas intialization has not been implemented yet ...');
    end
    % ============================================================== %
end

if ( Global_flag & Local_flag )
    
    % ============================================================== %
    % global part    
    % ============================================================== %
    if ( Kmeans_flag_Global )
             
        % detect PVs                          
        if ( classnumber == 5 )
            % run segmentation scripts
            if ( MRF_flag_Global )
                if ( GMM_PVs_flag )
                    Global_Kmeans_EM_MRF_5classes_GmmPVs;
                else
                    Global_Kmeans_EM_MRF_5classes;
                end
            else
                if ( GMM_PVs_flag )
                    Global_Kmeans_EM_5classes_PVs;
                else
                    Global_Kmeans_EM_5classes;
                end
            end
        end
        
        if ( classnumber == 4 )
            % run segmentation scripts
            if ( MRF_flag_Global )
                if ( GMM_PVs_flag )
                    Global_Kmeans_EM_MRF_4classes_GmmPVs;
                else
                    Global_Kmeans_EM_MRF_4classes;
                end
%                 Global_Kmeans_EM_MRF_4classes;
            else
                if ( GMM_PVs_flag )
                    Global_Kmeans_EM_4classes_PVs;
                else
                    Global_Kmeans_EM_4classes;
                end
            end
        end
    end
    % ============================================================== %
    if ( Atlas_flag_Global )
        disp('Atlas intialization has not been implemented yet ...');
    end
    % ============================================================== %
    
    % ============================================================== %
    % Local part    
    % ============================================================== %
    if ( Kmeans_flag_Local )
             
        % detect PVs                          
        if ( classnumber == 5 )
            % run segmentation scripts
            if ( MRF_flag_Local )
                if ( GMM_PVs_flag_Local )
                    Local_Kmeans_EM_MRF_5classes_GmmPVs; % all five classes
                else
                    Local_Kmeans_EM_MRF_5classes;
                end
            else
                if ( GMM_PVs_flag_Local )
                    Local_Kmeans_EM_5classes_PVs; % only gm + wm1 + wm2
                else
                    Local_Kmeans_EM_5classes;
                end
            end
        end
        
        if ( classnumber == 4 )
            % run segmentation scripts
            if ( MRF_flag_Local )
                if ( GMM_PVs_flag_Local )
                    Local_Kmeans_EM_MRF_4classes_GmmPVs;  % all four classes
                else
                    Local_Kmeans_EM_MRF_4classes;
                end
            else
                if ( GMM_PVs_flag_Local )
                    Local_Kmeans_EM_4classes_PVs; % only gm + wm
                else
                    Local_Kmeans_EM_4classes;  % only gm + wm
                end
            end
        end
    end
    % ============================================================== %
    if ( Atlas_flag_Local )
        disp('Atlas intialization has not been implemented yet ...');
    end
    % ============================================================== %
end