
AllinOne_ResetAll

disp('----------------------------------------------------');
disp('AllinOne_ComputeCortexStatistics_WithStem');

minThickness = 1.5;
datacolor = [1 1 1];
opacity = 1;
%==========================================================
% step 1: mean curvature

BrainDir = 'cortex_reconstruction';
brainmask_nostem_roi = 'brainmask_nostem_roi.hdr';

[subdirs, num] = FindAllDirectory(home)

for i=1:num
    dirName = fullfile(home, subdirs{i}, BrainDir)

    cd(dirName);

    filename = 'Middle_GaussCurvature.hdr';
    if ( isempty(dir(filename))==0 )
        continue;
    end

    %------------------------------------------------------------------
    filename = 'N3Brain_roi.hdr';
    p = dir(filename);
    if ( isempty(p) == 1 )
        continue;
    end
    [imagedata, header] = LoadAnalyze(filename, 'Grey');

    filename = 'Internal_levelset_Result.hdr';
    p = dir(filename);
    if ( isempty(p) == 0 )    
        [Internal, header] = LoadAnalyze(filename, 'Real');
    else
        Internal = [];
    end

    filename = 'External_levelset_Result_TH.hdr';
    p = dir(filename);
    if ( isempty(p) == 0 )
        [External, header] = LoadAnalyze(filename, 'Real');
    else
        External = [];
    end

%     filename = brainmask_nostem_roi;
%     p = dir(filename);
%     if ( isempty(p) == 0 )
%         [brainmask_nostem, header] = LoadAnalyze(filename, 'Grey');
%     else
%         brainmask_nostem = [];
%     end

    xsize = header.xsize;
    ysize = header.ysize;
    zsize = header.zsize;

    xvoxelsize = header.xvoxelsize;
    yvoxelsize = header.yvoxelsize;
    zvoxelsize = header.zvoxelsize;

    % perform the level set propagation
    %--------------------------------------------------------------------------
    accuracy = 'medium';
    tMax_ReIntialize = 100;
    errorMax = 0.05;
    %--------------------------------------------------------------------------
    CortexThickness_SignedDistanceFunction

    filename = 'Middle_levelset_Result.hdr';
    [Middle, header] = LoadAnalyze(filename, 'Real');

    % perform the level set propagation
    %--------------------------------------------------------------------------
    MeanCurvature_Cortex
end

%==========================================================
% step 2: meanTH, Std_TH, mean_InternalCurvature,std_InternalCurvature, surface area, cortex volume, IPR

cd(home)
filename = 'All_MeanThickness_MeanCurvature.mat';
brainmask_nostem_roi = 'brainmask_nostem_roi.hdr';
if ( isempty(dir(filename)) )

    CortexDir = 'cortex_reconstruction';

    level = 0;
%     minThickness = 0.5;

    N = 500;
    % name, meanTH, Std_TH, mean_InternalCurvature, std_InternalCurvature, surface area, cortex volume, IPR
    All_MeanThickness_MeanCurvature = cell(N, 8); 

    index = 0;

    [subdirs, num] = FindAllDirectory(home)

    for i=1:num
        dirName = fullfile(home, subdirs{i}, CortexDir)

        cd(dirName);
        
        filename = 'Internal_Cortex_Statistics.mat';
        if ( isempty(dir(filename))==1 )
            %------------------------------------------------------------------
            filename = 'N3Brain_roi.hdr';
            p = dir(filename);
            if ( isempty(p) )
                continue;
            end
            [imagedata, header] = LoadAnalyze(filename, 'Grey');

            filename = 'Internal_levelset_Result.hdr';
            p = dir(filename);
            if ( isempty(p) )
                continue;
            end
            [Internal, header] = LoadAnalyze(filename, 'Real');

            filename = 'External_levelset_Result_TH.hdr';
            p = dir(filename);
            if ( isempty(p) )
                continue;
            end
            [External, header] = LoadAnalyze(filename, 'Real');

            % -----------------------------------------------------------------
            filename = brainmask_nostem_roi;
            p = dir(filename);
            if ( isempty(p) == 0 )
                [brainmask_nostem, header] = LoadAnalyze(filename, 'Grey');
            else
                brainmask_nostem = [];
            end
            % -----------------------------------------------------------------

            filename = 'Internal_Curvature.hdr';
            p = dir(filename);
            if ( isempty(p) )
                continue;
            end
            [InternalCurvature, header] = LoadAnalyze(filename, 'Real');

            filename = 'Cortex_Thickness.hdr';
            p = dir(filename);
            if ( isempty(p) )
                continue;
            end
            [TH, header] = LoadAnalyze(filename, 'Real');

            xsize = header.xsize;
            ysize = header.ysize;
            zsize = header.zsize;

            xvoxelsize = header.xvoxelsize;
            yvoxelsize = header.yvoxelsize;
            zvoxelsize = header.zvoxelsize;

            voxelsize = (xvoxelsize+yvoxelsize+zvoxelsize)/3;
            InternalCurvature = abs(InternalCurvature);

            MeanThickness_Curvature_Cortex;
        else
            load(filename);
        end

        index = index + 1;
        All_MeanThickness_MeanCurvature{index, 1} = subdirs{i};
        All_MeanThickness_MeanCurvature{index, 2} = mean_TH;
        All_MeanThickness_MeanCurvature{index, 3} = std_TH;
        All_MeanThickness_MeanCurvature{index, 4} = mean_InternalCurvature;
        All_MeanThickness_MeanCurvature{index, 5} = std_InternalCurvature;

        All_MeanThickness_MeanCurvature{index, 6} = surface_area;
        All_MeanThickness_MeanCurvature{index, 7} = cortex_volume;        
        All_MeanThickness_MeanCurvature{index, 8} = IPR;        

    end

    All_MeanThickness_MeanCurvature = All_MeanThickness_MeanCurvature(1:index, :);
    cd(home);
    save All_MeanThickness_MeanCurvature All_MeanThickness_MeanCurvature

    All_MeanThickness_MeanCurvature

end
%==========================================================

% % step 3: mask the four lobes
% AllinOne_PropagatingFourLobes

%==========================================================

% % step 4: second order statistics
% CortexDir = 'cortex_reconstruction';
% brainmask_nostem_roi = 'brainmask_nostem_roi.hdr';
% 
% level = 0;
% % minThickness = 0.5;
% 
% N = 500;
% % name, ICI, GLN, ECI, MLN, CR, ICI_vtk, GLN_vtk, ECI_vtk, MLN_vtk, volume_CH
% All_SecondOrder_Statistics = cell(N, 11); 
% 
% % name, ICI, GLN, ECI, MLN, surface area
% Frontal_SecondOrder_Statistics = cell(N, 6);
% Occipital_SecondOrder_Statistics = cell(N, 6);
% Parietal_SecondOrder_Statistics = cell(N, 6);
% Temporal_SecondOrder_Statistics = cell(N, 6);
% 
% index = 0;
% 
% [subdirs, num] = FindAllDirectory(home)
% 
% for i=1:num
%     dirName = fullfile(home, subdirs{i}, CortexDir)
% 
%     cd(dirName);
%     
%     filename = 'Internal_SecondOrder_Cortex_Statistics.mat';
% %     if ( isempty(dir(filename)) )
%         %------------------------------------------------------------------
%         filename = 'N3Brain_roi.hdr';
%         p = dir(filename);
%         if ( isempty(p) )
%             continue;
%         end
%         [imagedata, header] = LoadAnalyze(filename, 'Grey');
% 
%         filename = 'Internal_levelset_Result.hdr';
%         p = dir(filename);
%         if ( isempty(p) )
%             continue;
%         end
%         [Internal, header] = LoadAnalyze(filename, 'Real');
% 
%         filename = 'External_levelset_Result_TH.hdr';
%         p = dir(filename);
%         if ( isempty(p) )
%             continue;
%         end
%         [External, header] = LoadAnalyze(filename, 'Real');
% 
%         % -----------------------------------------------------------------
%         filename = brainmask_nostem_roi;
%         p = dir(filename);
%         if ( isempty(p) == 0 )
%             [brainmask_nostem, header] = LoadAnalyze(filename, 'Grey');
%         else
%             brainmask_nostem = [];
%         end
%         % -----------------------------------------------------------------
% 
%         filename = 'Internal_Curvature.hdr';
%         p = dir(filename);
%         if ( isempty(p) )
%             continue;
%         end
%         [InternalCurvature, header] = LoadAnalyze(filename, 'Real');
% 
%         filename = 'Internal_GaussCurvature.hdr';
%         p = dir(filename);
%         if ( isempty(p) )
%             continue;
%         end
%         [InternalGaussCurvature, header] = LoadAnalyze(filename, 'Real');
% 
%         xsize = header.xsize;
%         ysize = header.ysize;
%         zsize = header.zsize;
% 
%         xvoxelsize = header.xvoxelsize;
%         yvoxelsize = header.yvoxelsize;
%         zvoxelsize = header.zvoxelsize;
% 
%         voxelsize = (xvoxelsize+yvoxelsize+zvoxelsize)/3;
% 
%         %MeanThickness_Curvature_Cortex;
%         SecondOrder_Statistics
% %     else
% %         load(filename);
% %     end
%     
%     index = index + 1;
%     All_SecondOrder_Statistics{index, 1} = subdirs{i};
%     All_SecondOrder_Statistics{index, 2} = ICI;
%     All_SecondOrder_Statistics{index, 3} = GLN;
%     All_SecondOrder_Statistics{index, 4} = ECI;
%     All_SecondOrder_Statistics{index, 5} = MLN;
%     All_SecondOrder_Statistics{index, 6} = CR;
% 
%     All_SecondOrder_Statistics{index, 7} = ICI_vtk;
%     All_SecondOrder_Statistics{index, 8} = GLN_vtk;
%     All_SecondOrder_Statistics{index, 9} = ECI_vtk;
%     All_SecondOrder_Statistics{index, 10} = MLN_vtk;
%     
%     All_SecondOrder_Statistics{index, 11} = volume_CH; % volume of QHull in ml
% 
%     Frontal_SecondOrder_Statistics{index, 1} = subdirs{i};
%     Frontal_SecondOrder_Statistics{index, 2} = ICI_frontal;
%     Frontal_SecondOrder_Statistics{index, 3} = GLN_frontal;
%     Frontal_SecondOrder_Statistics{index, 4} = ECI_frontal;
%     Frontal_SecondOrder_Statistics{index, 5} = MLN_frontal;
%     Frontal_SecondOrder_Statistics{index, 6} = SA_frontal;
% 
%     Occipital_SecondOrder_Statistics{index, 1} = subdirs{i};
%     Occipital_SecondOrder_Statistics{index, 2} = ICI_occipital;
%     Occipital_SecondOrder_Statistics{index, 3} = GLN_occipital;
%     Occipital_SecondOrder_Statistics{index, 4} = ECI_occipital;
%     Occipital_SecondOrder_Statistics{index, 5} = MLN_occipital;
%     Occipital_SecondOrder_Statistics{index, 6} = SA_occipital;
% 
%     Parietal_SecondOrder_Statistics{index, 1} = subdirs{i};
%     Parietal_SecondOrder_Statistics{index, 2} = ICI_parietal;
%     Parietal_SecondOrder_Statistics{index, 3} = GLN_parietal;
%     Parietal_SecondOrder_Statistics{index, 4} = ECI_parietal;
%     Parietal_SecondOrder_Statistics{index, 5} = MLN_parietal;
%     Parietal_SecondOrder_Statistics{index, 6} = SA_parietal;
% 
%     Temporal_SecondOrder_Statistics{index, 1} = subdirs{i};
%     Temporal_SecondOrder_Statistics{index, 2} = ICI_temporal;
%     Temporal_SecondOrder_Statistics{index, 3} = GLN_temporal;
%     Temporal_SecondOrder_Statistics{index, 4} = ECI_temporal;
%     Temporal_SecondOrder_Statistics{index, 5} = MLN_temporal;
%     Temporal_SecondOrder_Statistics{index, 6} = SA_temporal;
% 
% end
% All_SecondOrder_Statistics = All_SecondOrder_Statistics(1:index, :);
% Frontal_SecondOrder_Statistics = Frontal_SecondOrder_Statistics(1:index, :);
% Occipital_SecondOrder_Statistics = Occipital_SecondOrder_Statistics(1:index, :);
% Parietal_SecondOrder_Statistics = Parietal_SecondOrder_Statistics(1:index, :);
% Temporal_SecondOrder_Statistics = Temporal_SecondOrder_Statistics(1:index, :);
% 
% cd(home);
% save All_SecondOrder_Statistics All_SecondOrder_Statistics Frontal_SecondOrder_Statistics Occipital_SecondOrder_Statistics ...
%     Parietal_SecondOrder_Statistics Temporal_SecondOrder_Statistics
% 
% All_SecondOrder_Statistics

disp('----------------------------------------------------');

%==========================================================
% step 5: brain volume

cd(home)
filename = 'BrainVolume_results.mat';

% if ( isempty(dir(filename)) )

    brainmaskDir = 'brainMask';
    brainmask_noStem = 'brainmask_nostem.hdr';

    CortexDir = 'cortex_reconstruction';
    corpus_callosumfile = 'corpus_callosum.hdr';
    basal_gangliafile = 'basal_ganglia.hdr';

    N = 500;
    BrainVolume_results = cell(N, 5); % name, brain volume ml, gray matter volume, white matter volume, csf volume

    tt = 0;

    [subdirs, num] = FindAllDirectory(home)

    for i=1:num

        dirName = fullfile(home, subdirs{i})

        cd(dirName);

        disp('brain volume ...');
        
        filename = 'brain_volume_statistics.mat';
        p = dir(filename);
        %if ( isempty(p)==1 )
            %------------------------------------------------------------------    
            filename = 'segResult_4classes.hdr';
            p = dir(filename);
            if ( isempty(p) )
                continue;
            end
            [segResult, header] = LoadAnalyze(filename, 'Grey');        
            %------------------------------------------------------------------
            brainmaskName = fullfile(home, subdirs{i}, brainmaskDir, brainmask_noStem);
            if ( isempty(dir(brainmaskName))==1 )
                continue;
            end

            [brainmask, header] = LoadAnalyze(brainmaskName, 'Grey');  
            %------------------------------------------------------------------
            corpus_callosumName = fullfile(home, subdirs{i}, CortexDir, corpus_callosumfile);
            if ( isempty(dir(corpus_callosumName))==1 )
                continue;
            end

            [corpus_callosum, header] = LoadAnalyze(corpus_callosumName, 'Grey');            
            %------------------------------------------------------------------
            basal_gangliaName = fullfile(home, subdirs{i}, CortexDir, basal_gangliafile);
            if ( isempty(dir(basal_gangliaName))==1 )
                continue;
            end

            [basal_ganglia, header] = LoadAnalyze(basal_gangliaName, 'Grey');
            %------------------------------------------------------------------
            segResult(find(brainmask==0)) = 0;

            num_basal_ganglia = length(find(basal_ganglia>0));
            num_corpus_callosum = length(find(corpus_callosum>0));
            num_GrayMatter = length(find(segResult==2));
            num_WhiteMatter = length(find(segResult==3));
            num_CSF = length(find(segResult==1));

            VL_ratio = header.xvoxelsize*header.yvoxelsize*header.zvoxelsize* 10^-3;

            brainvolume = (num_GrayMatter+num_WhiteMatter+num_basal_ganglia)*VL_ratio; % mL
            graymattervolume = (num_GrayMatter-num_corpus_callosum)*VL_ratio; % mL
            whitemattervolume = (num_WhiteMatter+num_corpus_callosum)*VL_ratio; % mL
            csfvolume = num_CSF*VL_ratio; % mL

            save brain_volume_statistics brainvolume graymattervolume whitemattervolume csfvolume
            %------------------------------------------------------------------
%         else
%             load(filename);
%         end
        tt = tt+1;
        BrainVolume_results{tt,1} = subdirs{i};
        BrainVolume_results{tt,2} = brainvolume;
        BrainVolume_results{tt,3} = graymattervolume;
        BrainVolume_results{tt,4} = whitemattervolume;
        BrainVolume_results{tt,5} = csfvolume;
    end
    BrainVolume_results = BrainVolume_results(1:tt, :);
    cd(home);
    save BrainVolume_results BrainVolume_results
% end

%==========================================================
% step 6: volumes for four lobes

cd(home)
filename = 'FourLobes_volume_results.mat';

if ( isempty(dir(filename)) )

    brainmaskDir = 'brainMask';
    brainmask_noStem = 'brainmask_nostem.hdr';

    CortexDir = 'cortex_reconstruction';
    corpus_callosumfile = 'corpus_callosum.hdr';
    frontal_template = 'frontal_template.hdr';
    parietal_template = 'parietal_template.hdr';
    occipital_template = 'occipital_template.hdr';
    temporal_template = 'temporal_template.hdr';

    N = 500;
    FourLobes_results = cell(N, 9); 
    % name, frontal volume ml, parietal volume, occipital volume, temporal volume, frontal WM volume ml, parietal WM volume, occipital WM volume, temporal WM volume

    tt = 0;

    [subdirs, num] = FindAllDirectory(home)

    for i=1:num

        dirName = fullfile(home, subdirs{i})

        cd(dirName);

        disp('brain volume ...');
        
        filename = 'FourLobes_volume_statistics.mat';
        p = dir(filename);
%         if ( isempty(p)==1 )
            %------------------------------------------------------------------    
            filename = 'segResult_4classes.hdr';
            p = dir(filename);
            if ( isempty(p) )
                continue;
            end
            [segResult, header] = LoadAnalyze(filename, 'Grey');        
            %------------------------------------------------------------------
            brainmaskName = fullfile(home, subdirs{i}, brainmaskDir, brainmask_noStem);
            if ( isempty(dir(brainmaskName))==1 )
                continue;
            end

            [brainmask, header] = LoadAnalyze(brainmaskName, 'Grey');  
            %------------------------------------------------------------------
            corpus_callosumName = fullfile(home, subdirs{i}, CortexDir, corpus_callosumfile);
            if ( isempty(dir(corpus_callosumName))==1 )
                continue;
            end

            [corpus_callosum, header] = LoadAnalyze(corpus_callosumName, 'Grey');            
            %------------------------------------------------------------------
            frontalName = fullfile(home, subdirs{i}, CortexDir, frontal_template);
            if ( isempty(dir(frontalName))==1 )
                continue;
            end

            [frontal_lobe, header] = LoadAnalyze(frontalName, 'Grey');            
            %------------------------------------------------------------------
            parietalName = fullfile(home, subdirs{i}, CortexDir, parietal_template);
            if ( isempty(dir(parietalName))==1 )
                continue;
            end

            [parietal_lobe, header] = LoadAnalyze(parietalName, 'Grey');            
            %------------------------------------------------------------------
            occipitalName = fullfile(home, subdirs{i}, CortexDir, occipital_template);
            if ( isempty(dir(occipitalName))==1 )
                continue;
            end

            [occipital_lobe, header] = LoadAnalyze(occipitalName, 'Grey');            
            %------------------------------------------------------------------
            temporalName = fullfile(home, subdirs{i}, CortexDir, temporal_template);
            if ( isempty(dir(temporalName))==1 )
                continue;
            end

            [temporal_lobe, header] = LoadAnalyze(temporalName, 'Grey');            
            %------------------------------------------------------------------

            segResult(find(brainmask==0)) = 0;
            
            GM_seg = segResult;
            GM_seg(find(GM_seg~=2))=0;
            GM_seg(find(corpus_callosum>0))=0;

            WM_seg = segResult;
            WM_seg(find(WM_seg~=3))=0;
            WM_seg(find(corpus_callosum>0))=1;
            
            num_frontal = length(find( (frontal_lobe>0) & (GM_seg>0) ));
            num_parietal = length(find( (parietal_lobe>0) & (GM_seg>0) ));
            num_occipital = length(find( (occipital_lobe>0) & (GM_seg>0) ));
            num_temporal = length(find( (temporal_lobe>0) & (GM_seg>0) ));

            VL_ratio = header.xvoxelsize*header.yvoxelsize*header.zvoxelsize* 10^-3;

            frontal_volume = num_frontal*VL_ratio; % mL
            parietal_volume = num_parietal*VL_ratio; % mL
            occipital_volume = num_occipital*VL_ratio; % mL
            temporal_volume = num_temporal*VL_ratio; % mL

            num_frontal = length(find( (frontal_lobe>0) & (WM_seg>0) ));
            num_parietal = length(find( (parietal_lobe>0) & (WM_seg>0) ));
            num_occipital = length(find( (occipital_lobe>0) & (WM_seg>0) ));
            num_temporal = length(find( (temporal_lobe>0) & (WM_seg>0) ));

            VL_ratio = header.xvoxelsize*header.yvoxelsize*header.zvoxelsize* 10^-3;

            frontal_volume_WM = num_frontal*VL_ratio; % mL
            parietal_volume_WM = num_parietal*VL_ratio; % mL
            occipital_volume_WM = num_occipital*VL_ratio; % mL
            temporal_volume_WM = num_temporal*VL_ratio; % mL

            save FourLobes_volume_statistics frontal_volume parietal_volume occipital_volume temporal_volume ...
                frontal_volume_WM parietal_volume_WM occipital_volume_WM temporal_volume_WM
            %------------------------------------------------------------------
%         else
%             load(filename);
%         end
        tt = tt+1;
        FourLobes_results{tt,1} = subdirs{i};
        FourLobes_results{tt,2} = frontal_volume;
        FourLobes_results{tt,3} = parietal_volume;
        FourLobes_results{tt,4} = occipital_volume;
        FourLobes_results{tt,5} = temporal_volume;
        
        FourLobes_results{tt,6} = frontal_volume_WM;
        FourLobes_results{tt,7} = parietal_volume_WM;
        FourLobes_results{tt,8} = occipital_volume_WM;
        FourLobes_results{tt,9} = temporal_volume_WM;

    end
    FourLobes_volume_results = FourLobes_results(1:tt, :);
    cd(home);
    save FourLobes_volume_results FourLobes_volume_results
end

%==========================================================

% % step 7: four lobes statistics
% CortexDir = 'cortex_reconstruction';
% brainmask_nostem_roi = 'brainmask_nostem_roi.hdr';
% 
% N = 500;
% 
% % name, Surface area, cortical volume, mean Thickness, CR (convexity ratio), QH_volumes
% Frontal_Statistics = cell(N, 6);
% Occipital_Statistics = cell(N, 6);
% Parietal_Statistics = cell(N, 6);
% Temporal_Statistics = cell(N, 6);
% 
% index = 0;
% 
% [subdirs, num] = FindAllDirectory(home)
% 
% for i=1:num
%     dirName = fullfile(home, subdirs{i}, CortexDir)
% 
%     cd(dirName);
%     filename = 'FourLobes_Statistics';
% %     if ( isempty(dir(filename)) )
%         %------------------------------------------------------------------
%         filename = 'N3Brain_roi.hdr';
%         p = dir(filename);
%         if ( isempty(p) )
%             continue;
%         end
%         [imagedata, header] = LoadAnalyze(filename, 'Grey');
% 
%         xsize = header.xsize;
%         ysize = header.ysize;
%         zsize = header.zsize;
% 
%         xvoxelsize = header.xvoxelsize;
%         yvoxelsize = header.yvoxelsize;
%         zvoxelsize = header.zvoxelsize;
% 
%         voxelsize = (xvoxelsize+yvoxelsize+zvoxelsize)/3;
% 
%         FourLobes_Statistics
% %     else
% %         load(filename);
% %     end
%     
%     index = index + 1;
%     
%     % name, Surface area, cortical volume, mean Thickness, CR (convexity ratio)
% 
%     Frontal_Statistics{index, 1} = subdirs{i};
%     Frontal_Statistics{index, 2} = SA_frontal;
%     Frontal_Statistics{index, 3} = CV_frontal;
%     Frontal_Statistics{index, 4} = MeanTH_frontal;
%     Frontal_Statistics{index, 5} = CR_frontal;
%     Frontal_Statistics{index, 6} = frontal_QH_volume;
% 
%     Occipital_Statistics{index, 1} = subdirs{i};
%     Occipital_Statistics{index, 2} = SA_occipital;
%     Occipital_Statistics{index, 3} = CV_occipital;
%     Occipital_Statistics{index, 4} = MeanTH_occipital;
%     Occipital_Statistics{index, 5} = CR_occipital;
%     Occipital_Statistics{index, 6} = occipital_QH_volume;
% 
%     Parietal_Statistics{index, 1} = subdirs{i};
%     Parietal_Statistics{index, 2} = SA_parietal;
%     Parietal_Statistics{index, 3} = CV_parietal;
%     Parietal_Statistics{index, 4} = MeanTH_parietal;
%     Parietal_Statistics{index, 5} = CR_parietal;
%     Parietal_Statistics{index, 6} = parietal_QH_volume;
% 
%     Temporal_Statistics{index, 1} = subdirs{i};
%     Temporal_Statistics{index, 2} = SA_temporal;
%     Temporal_Statistics{index, 3} = CV_temporal;
%     Temporal_Statistics{index, 4} = MeanTH_temporal;
%     Temporal_Statistics{index, 5} = CR_temporal;
%     Temporal_Statistics{index, 6} = temporal_QH_volume;
%    
% end
% Frontal_Statistics = Frontal_Statistics(1:index, :);
% Occipital_Statistics = Occipital_Statistics(1:index, :);
% Parietal_Statistics = Parietal_Statistics(1:index, :);
% Temporal_Statistics = Temporal_Statistics(1:index, :);

% cd(home);
% save All_FourLobes_Statistics Frontal_Statistics Occipital_Statistics ...
%     Parietal_Statistics Temporal_Statistics

disp('----------------------------------------------------');