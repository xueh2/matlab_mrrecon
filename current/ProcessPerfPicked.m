
clear all
closeall

data_dir = '/fastdata/data/perf//'

data_dir = '/fastdata/data/highres/'


findAndMoveMeasDat(data_dir)

UTCases = set_up_UT_cases_Perfusion_PickedCase

[names, num] = FindSubDirs(data_dir)


for i=1:num
    try
        UTCases{1, 1} = data_dir;
        UTCases{1, 2} = names{i}

        UTCases{1, 4} = 'GT_QPerf_AI_STCNNT_OFFLINE.xml'
        UTCases{1, 5} = 'res_GT_QPerf_AI_STCNNT_OFFLINE' 

        % UTCases{1, 4} = 'GT_QPerf_AI_STCNNT_istore.xml'
        % UTCases{1, 5} = 'res_GT_QPerf_AI_STCNNT_istore' 
        % 
        % UTCases{1, 4} = 'GT_QPerf_AI.xml'
        % UTCases{1, 5} = 'res_GT_QPerf_AI' 


        if strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0
            UTCases{1, 4} = 'default_measurement_dependencies_Noise_CoilSen_SCC.xml'
        end

        performUTValidation(UTCases, 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, UTCases{1, 4}, [], 'IsmrmrdParameterMap_Siemens_Perfusion_NX.xml', '/home/xueh/Debug/DebugOutput')
    
        %performUTValidation(UTCases, 0, 0, 'GCRSANDBOX455.redmond.corp.microsoft.com', '9002', 1, 1, 0, 0, 0, UTCases{1, 4}, [], [], '/home/xueh/Debug/DebugOutput')

        %performUTValidation(UTCases, 0, 0, 'GCRAZGDL3156.westus2.cloudapp.azure.com', '9002', 1, 1, 0, 0, 0, UTCases{1, 4}, [], [], '/home/xueh/Debug/DebugOutput')

        debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput')

    catch
    end
end

mat_dir = '/data/raw_data/highres_mat/'
mkdir(mat_dir)

for i=1:num
    try
        UTCases{1, 1} = data_dir;
        UTCases{1, 2} = names{i}

        UTCases{1, 4} = 'GT_QPerf_AI_STCNNT_OFFLINE.xml'
        UTCases{1, 5} = 'res_GT_QPerf_AI_STCNNT_OFFLINE' 

        if isempty(strfind(names{i}, 'ISMRMRD_Noise_dependency'))

            try
                debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput')
                [data, gmap, res] = prepare_for_stcnnt_inference_perf(debug_dir); 
                h1 = figure; imagescn(cat(4, data, res, data-res), [], [3 size(data, 4)], [8], 3);
                saveas(h1, fullfile(data_dir, [names{i} '_raw_ai.fig']));
                fmaps = load_maps_for_perf(debug_dir);
                f1 = figure('Name', names{i}); imagescn(fmaps, [0 8], [1 size(fmaps, 3)], 12); PerfColorMap
                saveas(f1, fullfile(data_dir, [names{i} '_fmaps_GT_QPerf_AI_STCN NT_OFFLINE_with_pre_mapping.fig']));

                % Q_e = analyze75read(fullfile(debug_dir, 'Q_e.hdr'));
                % gd_0 = analyze75read(fullfile(debug_dir, 'CASignal_Perf_PSIR_0.hdr'));
                % gd_1 = analyze75read(fullfile(debug_dir, 'CASignal_Perf_PSIR_1.hdr'));
                % gd_2 = analyze75read(fullfile(debug_dir, 'CASignal_Perf_PSIR_2.hdr'));
                % 
                % save(fullfile(mat_dir, [names{i} '.mat']), 'data', 'gmap', 'res', 'fmaps', 'Q_e', 'gd_0', 'gd_1', 'gd_2');
            catch
            end
        end
    
        closeall
    catch
    end
end

for i=1:num
    if strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0
        continue
    end
    disp([num2str(i) ' --- ' names{i}])
    try
        open(fullfile(data_dir, [names{i} '_fmaps_GT_QPerf_AI_STCNNT_OFFLINE_with_pre_mapping.fig']));
        MR_toolbox        
        pause
        closeall
    catch
    end
end

res_dir = '/data/raw_data/highres/Perfusion_AIF_Q_mapping_66016_001246866_001246875_267_20250127-123612/res_GT_QPerf_AI_STCNNT_OFFLINE/DebugOutput//model_res_73'

res_dir = '/data/raw_data/highres/Perfusion_AIF_Q_mapping_66016_111101640_111101649_388_20250306-141411/res_GT_QPerf_AI_STCNNT_OFFLINE/DebugOutput/model_res_73'

[input, res, gmap] = load_results_stcnnt_inference_perf(res_dir, 1);

%% check perf debug

cd /data/raw_data/perf_debug/scc/Perfusion_AIF_TwoEchoes_Interleaved_R2_66016_31820219_31820228_802_20241119-150332/res_GT_QPerf_AI_STCNNT_OFFLINE/DebugOutput/

a = analyze75read('compute_final_btex__final_F_job_0.hdr');
figure; imagescn(cat(3, a, a/(1-0.42)), [0 8], [], [16]); PerfColorMap;

b = analyze75read('compute_final_btex__final_delay_job_0.hdr');
figure; imagescn(b, [0 3], [], [16]); PerfColorMap;

Q_e = analyze75read('compute_final_btex__fitted_Q_e_job_0.hdr');
size(Q_e)

Gd = analyze75read('compute_cost_Gd_vs_delayed_Q_e__Gd_job_0.hdr');
size(Gd)

Gd2 = analyze75read('PerfFlowMapping_Job_0_perf_moco_upsampled.hdr');
size(Gd2)

figure; imagescn(cat(4, Gd2, Gd, Q_e), [0 1], [], [16], 3)

% -------------------------------------
%% spine

data_dir = '/data/raw_data/20250122_SPINE_DATFILES//'

data_dir = '/data/raw_data/noncardiac/spine/'

data_dir = '/data/raw_data/noncardiac_20250805/knee/'

data_dir = '/data/raw_data/noncardiac/phantom/'

data_dir = '/fastdata/noncardiac/phantom/20250815//'

data_dir = '/home/gtuser/data/noncardiac/'

findAndMoveMeasDat(data_dir)

UTCases = set_up_UT_cases_knee;

[names, num] = FindSubDirs(data_dir)

for i=1:num
    UTCases{1, 1} = data_dir;
    UTCases{1, 2} = names{i}

    if strfind(names{i}, 'Spine') > 0
        UTCases{1, 4} = 'GTPrep_2DT_STCNNT_Spine.xml'
        UTCases{1, 5} = 'res_GTPrep_2DT_STCNNT_Spine'
    end

    if strfind(names{i}, 'Data_3DT') > 0
        UTCases{1, 4} = 'Generic_Cartesian_3D_Grappa_STCNNT.xml'
        UTCases{1, 5} = 'res_Generic_Cartesian_3D_Grappa_STCNNT'
    end

    if strfind(names{i}, 'Data_2D') > 0
        UTCases{1, 4} = 'Generic_Cartesian_2D_Grappa_STCNNT.xml'
        UTCases{1, 5} = 'res_Generic_Cartesian_2D_Grappa_STCNNT'
    end

    h5name = fullfile(data_dir, names{i}, [names{i} '.h5']);
    if exist(h5name)
        dset = ismrmrd.Dataset(h5name);
        header = ismrmrd.xml.deserialize(dset.readxml());
        R = header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1
        dset.close()

        header.measurementInformation.protocolName
    end

    if any([(strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0) (strfind(names{i}, 'noise_') > 0)])
            UTCases{1, 4} = 'default_measurement_dependencies_Noise_CoilSen_SCC.xml'
    end

    UTCases{1, 4}
    performUTValidation(UTCases(1,:), 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], 'wip_070_fire_IsmrmrdParameterMap_Siemens.xml', '/home/gtuser/Debug/DebugOutput');
end

for i=1:num
    UTCases{1, 1} = data_dir;
    UTCases{1, 2} = names{i};

    if strfind(names{i}, 'Spine') > 0
        UTCases{1, 4} = 'GTPrep_2DT_STCNNT_Spine.xml'
        UTCases{1, 5} = 'res_GTPrep_2DT_STCNNT_Spine'
    end

    if strfind(names{i}, 'Data_3DT') > 0
        UTCases{1, 4} = 'Generic_Cartesian_3D_Grappa_STCNNT.xml'
        UTCases{1, 5} = 'res_Generic_Cartesian_3D_Grappa_STCNNT'
    end

    if strfind(names{i}, 'Data_2D') > 0
        UTCases{1, 4} = 'Generic_Cartesian_2D_Grappa_STCNNT.xml'
        UTCases{1, 5} = 'res_Generic_Cartesian_2D_Grappa_STCNNT'
    end

    h5name = fullfile(data_dir, names{i}, [names{i} '.h5']);
    if exist(h5name)
        dset = ismrmrd.Dataset(h5name);
        header = ismrmrd.xml.deserialize(dset.readxml());
        R = header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1
        dset.close()

        header.measurementInformation.protocolName

    end

    if any([(strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0) (strfind(names{i}, 'noise_') > 0)])
            continue
    end

    UTCases{1, 4}

    try
        debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
        [data, gmap, res] = prepare_for_stcnnt_inference_spine(debug_dir); 
        if ~isempty(res)
            SLC = size(data, 3)
            h = figure('Name', names{i}); imagescn(abs(cat(4, data, res)), [], [1 2], [8], 3);
            saveas(h, fullfile(data_dir, [names{i} '_' UTCases{1, 5} '.fig']));
        end
        save(fullfile(data_dir, [names{i} '.mat']), 'data', 'gmap', 'res', 'header');
    catch
    end
end

for i=1:num
    UTCases{1, 1} = data_dir;
    UTCases{1, 2} = names{i};

    if strfind(names{i}, 'Spine') > 0
        UTCases{1, 4} = 'GTPrep_2DT_STCNNT_Spine.xml';
        UTCases{1, 5} = 'res_GTPrep_2DT_STCNNT_Spine';
    end

    if strfind(names{i}, 'Data_3DT') > 0
        UTCases{1, 4} = 'Generic_Cartesian_3D_Grappa_STCNNT.xml';
        UTCases{1, 5} = 'res_Generic_Cartesian_3D_Grappa_STCNNT';
    end

    if strfind(names{i}, 'Data_2D') > 0
        UTCases{1, 4} = 'Generic_Cartesian_2D_Grappa_STCNNT.xml';
        UTCases{1, 5} = 'res_Generic_Cartesian_2D_Grappa_STCNNT';
    end

    h5name = fullfile(data_dir, names{i}, [names{i} '.h5']);
    if exist(h5name)
        dset = ismrmrd.Dataset(h5name);
        header = ismrmrd.xml.deserialize(dset.readxml());
        R = header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1;
        dset.close();

        header.measurementInformation.protocolName;

    end

    if any([(strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0) (strfind(names{i}, 'noise_') > 0)])
            continue
    end

    disp([names{i} ' - ' header.measurementInformation.protocolName ' - R = ' num2str(R) ' - matrix size ' num2str(header.encoding.reconSpace.matrixSize.x) ' ' num2str(header.encoding.reconSpace.matrixSize.y)])
end




res_dir = '/data/raw_data/noncardiac/spine/t1_tse_sag_p2_41837_59824231_59824236_450_20250804-172400/res_GTPrep_2DT_STCNNT_Spine/DebugOutput//model_res_62/'

[input, res, gmap] = load_results_stcnnt_inference_perf(res_dir, 1);


%% loc

data_dir = '/data/raw_data/debug/loc/'

data_dir = '/fastdata/phantom//'

findAndMoveMeasDat(data_dir)

UTCases = set_up_UT_cases_LGE;

[names, num] = FindSubDirs(data_dir)

names

for i=1:num
    UTCases{1, 1} = data_dir;
    UTCases{1, 2} = names{i}

    UTCases{1, 4} = 'Generic_Cartesian_Grappa_SNR.xml'
    UTCases{1, 5} = 'res_Generic_Cartesian_Grappa_SNR'     

    if strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0
        UTCases{1, 4} = 'default_measurement_dependencies_Noise_CoilSen_SCC.xml'
    end

    performUTValidation(UTCases(1,:), 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], 'IsmrmrdParameterMap_Siemens.xml', '/home/xueh/Debug/DebugOutput')
    names{i}

    debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
    [data, gmap, res] = prepare_for_stcnnt_inference_cine(debug_dir); 
    if ~isempty(res)
        SLC = size(gmap, 3)
        if SLC==1
            h = figure('Name', names{i}); imagescn(cat(4, data, res), [0 4*mean(abs(data(:)))], [SLC 2], [12], 3);
        else
            h = figure('Name', names{i}); imagescn(cat(4, data, res), [0 4*mean(abs(data(:)))], [2 SLC], [12], 3);
        end
    end
end

%% RetroCine

data_dir = '/data/raw_data/retro_cine/'

data_dir = '/data/raw_data/00035000/'

data_dir = '/data/raw_data/debug/cine/'

data_dir = '/fastdata/data/test'

data_dir = '/fastdata/data/debug/'

data_dir = '/data/raw_data/retrocine/00035000/'

findAndMoveMeasDat(data_dir)

UTCases = set_up_UT_cases_LGE;

[names, num] = FindSubDirs(data_dir)

names

for i=1:num
    UTCases{1, 1} = data_dir;
    UTCases{1, 2} = names{i}

    UTCases{1, 4} = 'GTPrep_2DT_RetroCine_STCNNT_offline.xml'
    UTCases{1, 5} = 'res_GTPrep_2DT_RetroCine_STCNNT_offline'     

    UTCases{1, 4} = 'GT_RetroCine_gmap_augmentation.xml'
    UTCases{1, 5} = 'res_GT_RetroCine_gmap_augmentation' 

    UTCases{1, 4} = 'GT_RetroCine.xml'
    UTCases{1, 5} = 'res_GT_RetroCine' 

    UTCases{1, 4} = 'GT_RetroCine_AI.xml'
    UTCases{1, 5} = 'res_GT_RetroCine_AI'

    UTCases{1, 4} = 'GTPrep_2DT_RetroCine_STCNNT_Dealiasing_SCC_GLS_Seg_AI_istore.xml'
    UTCases{1, 5} = 'res_GTPrep_2DT_RetroCine_STCNNT_Dealiasing_SCC_GLS_Seg_AI_istore'

    UTCases{1, 4} = 'GTPrep_2DT_RetroCine_SCC_GLS_Seg_AI_istore.xml'
    UTCases{1, 5} = 'res_GTPrep_2DT_RetroCine_SCC_GLS_Seg_AI_istore'

    if strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0
        UTCases{1, 4} = 'default_measurement_dependencies_Noise_CoilSen_SCC.xml'
    end

    performUTValidation(UTCases, 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput');
    names{i}

    % debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
    % [data, gmap, res] = prepare_for_stcnnt_inference_cine(debug_dir); 
    % if ~isempty(res)
    %     SLC = size(gmap, 3)
    %     if SLC==1
    %         h = figure('Name', names{i}); imagescn(cat(4, data, res), [0 4*mean(abs(data(:)))], [SLC 2], [12], 3);
    %     else
    %         h = figure('Name', names{i}); imagescn(cat(4, data, res), [0 4*mean(abs(data(:)))], [2 SLC], [12], 3);
    %     end
    % end
end

for i=1:num
    debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
    [data, gmap, res] = prepare_for_stcnnt_inference_cine(debug_dir); 
    if ~isempty(res)
        SLC = size(gmap, 3)
        if SLC==1
            h = figure('Name', names{i}); imagescn(cat(4, data, res), [], [SLC 2], [12], 3);
        else
            h = figure('Name', names{i}); imagescn(cat(4, data, res), [], [2 SLC], [12], 3);
        end
    end
end

for i=1:num
    [a, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 2100, 1);
    %[b, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 2101, 1);

    size(a)

    h = figure('name', names{i});
    SLC = size(a, 3)
    imagescn(a, [0 2.0*mean(a(:))], [1 SLC], [12], 4);

    saveas(h, fullfile(data_dir, [names{i} '.fig']), 'fig');
end

%% LGE

data_dir = '//data/raw_data/wb_psir/'

data_dir = '/data/raw_data/lge_flash//'

findAndMoveMeasDat(data_dir)

UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;

[names, num] = FindSubDirs(data_dir)

for i=1:num
    UTCases{1, 1} = data_dir;
    UTCases{1, 2} = names{i}

    for k=3:3
        if k == 1
            UTCases{1, 4} = 'GT_LGE.xml'
            UTCases{1, 5} = 'res_GT_LGE' 
        else 
            if k==2
                UTCases{1, 4} = 'GTPrep_WB_LGE_STCNNT.xml'
                UTCases{1, 5} = 'res_GTPrep_WB_LGE_STCNNT_old' 
            else 
                if k==3
                    UTCases{1, 4} = 'GTPrep_WB_LGE_STCNNT.xml'
                    UTCases{1, 5} = 'res_GTPrep_WB_LGE_STCNNT' 
                end
            end
        end
    
        if any([(strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0) (strfind(names{i}, 'noise_') > 0)])
                UTCases{1, 4} = 'default_measurement_dependencies.xml'
        end
    
        performUTValidation(UTCases, 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput');
        names{i}
    
        [a, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 111, 1);
    
        size(a)
    
        SLC = size(a, 3)
        for slc=1:SLC
            h = figure('name', names{i});
            imagescn(a(:,:,slc), [header(slc).window_center-header(slc).window_width/2 header(slc).window_center+header(slc).window_width/2], [], [8]);
       
            saveas(h, fullfile(data_dir, [names{i} '_' UTCases{1, 5} '_slc' num2str(slc-1) '_PSIR.fig']), 'fig');close(h)
        end
        

        try
            debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
            [data, gmap, res] = prepare_for_stcnnt_inference_LGE(debug_dir); 
            if ~isempty(res)
                SLC = size(data, 5)
                if SLC==1
                    h = figure('Name', names{i}); imagescn(cat(4, data, res), [], [SLC 4], [12], 3);
                else
                    h = figure('Name', names{i}); imagescn(cat(4, data, res), [], [4 SLC], [12], 3);
                end
                saveas(h, fullfile(data_dir, [names{i} '_' UTCases{1, 5} '.fig']));
                close(h)
            end
        catch
        end

        delete(fullfil('/tmp/gadgetron_data/*.h5'), 'f');
    end
end

res_dir = '/data/raw_data/lge_flash/LGE_MOCO_AVE_OnTheFly_66016_462800584_462800593_415_20250725-152130/res_GTPrep_WB_LGE_STCNNT/DebugOutput//model_res_70'
[input, res, gmap] = load_results_stcnnt_inference_perf(res_dir, 1);


for i=1:num
    if any([(strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0) (strfind(names{i}, 'noise_') > 0)])
        continue
    end

    try
        [a, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 111, 1);
    
        size(a)
    
        SLC = size(a, 3)
        for slc=1:SLC
            h = figure('name', names{i});
            imagescn(a(:,:,slc), [header(slc).window_center-header(slc).window_width/2 header(slc).window_center+header(slc).window_width/2], [], [8]);
       
            saveas(h, fullfile(data_dir, [names{i} '_' UTCases{1, 5} '_AI_slc' num2str(slc-1) '_PSIR.fig']), 'fig');
        end
        close(h)
    catch
        closeall;
    end
end

fig_dir = '/data/raw_data/wb_psir_fig/'
mkdir(fig_dir)

for i=1:num
    [a, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, 'res_GT_LGE'), 111, 1);
    b = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, 'res_GTPrep_WB_LGE_STCNNT_old'), 111);
    c = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, 'res_GTPrep_WB_LGE_STCNNT'), 111);

    size(a)

    SLC = size(a, 3)
    for slc=1:SLC
        h = figure('name', names{i});
        imagescn(cat(3, a(:,:,slc), b(:,:,slc), c(:,:,slc)), [header(slc).window_center-header(slc).window_width/2 header(slc).window_center+header(slc).window_width/2], [1 3], [16]);
    
        saveas(h, fullfile(fig_dir, [names{i} '_raw_AI_old_new_slc' num2str(slc-1) '_PSIR.fig']), 'fig');
        %close(h)
    end    
end


%% FatWater

data_dir = '/data/raw_data/FatWater/highres/'

UTCases = set_up_UT_cases_FatWater;

[names, num] = FindSubDirs(data_dir)

fat = {}
water = {}
for i=1:num
    UTCases{1, 1} = data_dir;
    UTCases{1, 2} = names{i}

    UTCases{1, 4} = 'GTPrep_2DT_FW_MOCO_AVE_STCNNT.xml'
    UTCases{1, 5} = 'res_GTPrep_2DT_FW_MOCO_AVE_STCNNT'

    if ~isempty(strfind(names{i}, 'PSIR'))
        UTCases{1, 4} = 'GTPrep_2DT_FW_MOCO_AVE_PSIR_STCNNT.xml'
        UTCases{1, 5} = 'res_GTPrep_2DT_FW_MOCO_AVE_PSIR_STCNNT'
    end

    % UTCases{1, 4} = 'GTPrep_2DT_FW_MOCO_AVE_Diego.xml'
    % UTCases{1, 5} = 'res_GTPrep_2DT_FW_MOCO_AVE_Diego'

    debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
    if exist(debug_dir)
        rmdir(debug_dir, 's');
    end

    performUTValidation(UTCases, 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
    names{i}
    
    try
        [data, gmap, res] = prepare_for_stcnnt_inference_FW(debug_dir);
        CON = size(res, 4)
        if CON == 1
            h = figure; imagescn(cat(4, data, res), [0 50], [CON 2], [18], 3);
        else
            h = figure; imagescn(cat(4, data, res), [0 50], [2 CON], [18], 3);
        end
    catch
    end

    fat{i} = readGTPlusExportData(fullfile(debug_dir, 'fat_CHA0_AVE0_SLC0_PHS0_REP0_0'));
    water{i} = readGTPlusExportData(fullfile(debug_dir, 'water_CHA0_AVE0_SLC0_PHS0_REP0_0'));
    h2 = figure; imagescn(cat(3, water{i}, fat{i}), [0 740], [], [18]);
    saveas(h2, fullfile(data_dir, [names{i} '_' UTCases{1, 5} '.fig']), 'fig')
end

fat = {}
water = {}
for i=1:num
    UTCases{1, 1} = data_dir;
    UTCases{1, 2} = names{i}
    UTCases{1, 4} = 'GTPrep_2DT_FW_MOCO_AVE_STCNNT.xml'
    UTCases{1, 5} = 'res_GTPrep_2DT_FW_MOCO_AVE_STCNNT' 

    % UTCases{1, 4} = 'GTPrep_2DT_FW_MOCO_AVE_Diego.xml'
    % UTCases{1, 5} = 'res_GTPrep_2DT_FW_MOCO_AVE_Diego'

    debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
    fat{i} = readGTPlusExportData(fullfile(debug_dir, 'fat_CHA0_AVE0_SLC0_PHS0_REP0_0'));
    water{i} = readGTPlusExportData(fullfile(debug_dir, 'water_CHA0_AVE0_SLC0_PHS0_REP0_0'));
end

f = Matlab_gt_resize_2D_image(double(abs(fat{1})), size(fat{2}, 1), size(fat{2}, 2), 5);
w = Matlab_gt_resize_2D_image(double(abs(water{1})), size(water{2}, 1), size(water{2}, 2), 5);

h = figure; imagescn(abs(cat(3, w, water{2}, water{3}, f, fat{2}, fat{3})), [0 740], [2 3], [18]);
saveas(h, fullfile(data_dir, ['R2_3_4_FW_' UTCases{1, 5} '.fig']), 'fig')

res_dir = '/data/raw_data/FatWater/highres/FatWater_MOCO_AVE_41837_2297927780_2297927785_389_20241216-173311/res_GTPrep_2DT_FW_MOCO_AVE_STCNNT/DebugOutput/model_res'

[input, res, gmap] = load_results_stcnnt_inference_perf(res_dir, 0);
CON = size(res, 4)
h = figure; imagescn(cat(4, data, res), [], [2 CON], [18], 3);

%% rtcine

data_dir = '/data/raw_data/rtcine/'

data_dir = '/data/raw_data/SNR_PAPER/SNR_PAPER/R5/'

data_dir = '/fastdata/data/realtime/'

data_dir = '/data/raw_data/SNR_PAPER/SNR_PAPER/R6/'

data_dir = '/home/gtuser/data/rtcine'


UTCases = set_up_UT_cases_RTCine;

findAndMoveMeasDat(data_dir)
[names, num] = FindSubDirs(data_dir)

for i=1:num
    UTCases{1, 1} = data_dir;
    UTCases{1, 2} = names{i}
   
    UTCases{1, 4} = 'Generic_RTCine_PInterp_Fil_ECG_STCNNT.xml'
    UTCases{1, 5} = 'res_Generic_RTCine_PInterp_Fil_ECG_STCNNT' 

    UTCases{1, 4} = 'Generic_RTCine_STCNNT.xml'
    UTCases{1, 5} = 'res_Generic_RTCine_STCNNT' 

    UTCases{1, 4} = 'GTPrep_2DT_RTCine.xml'
    UTCases{1, 5} = 'res_GTPrep_2DT_RTCine' 


     if any([(strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0) (strfind(names{i}, 'noise_') > 0)])
            UTCases{1, 4} = 'default_measurement_dependencies.xml'
     end

    performUTValidation(UTCases(1,:), 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
    names{i}

    % UTCases{1, 4} = 'Generic_RTCine.xml'
    % UTCases{1, 5} = 'res_Generic_RTCine' 

    try
        dset = ismrmrd.Dataset(fullfile(data_dir, names{i}, [names{i} '.h5']), 'dataset');
        hdr = ismrmrd.xml.deserialize(dset.readxml);
        dset.close()

        debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
        [data, gmap, res] = prepare_for_stcnnt_inference_cine(debug_dir); 
        h = figure('Name', [names{i} '_' hdr.measurementInformation.protocolName '_R' num2str(hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1) '_RO' num2str(hdr.encoding.reconSpace.matrixSize.x)]); 
        if size(data, 4) <= 2
            imagescn(cat(4, data, res), [0 3*abs(mean(data(:)))], [size(data, 4), 2], [12], 3);
        else
            imagescn(cat(4, data, res), [0 3*abs(mean(data(:)))], [2 size(data, 4)], [12], 3);
        end
        saveas(h, fullfile(data_dir, [names{i} '_' hdr.measurementInformation.protocolName '_R' num2str(hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1) '_RO' num2str(hdr.encoding.reconSpace.matrixSize.x) '.fig']));
    catch
    end
end

%% Fat Water

data_dir = '/data/raw_data/FatWater/00073013/'

[names, num] = FindSubDirs(data_dir)

UTCases = set_up_UT_cases_FatWater

for i=1:num
    UTCases{1, 1} = data_dir;
    UTCases{1, 2} = names{i}
    UTCases{1, 4} = 'GTPrep_2DT_FW_MOCO_AVE_PSIR_STCNNT.xml'
    UTCases{1, 5} = 'res_GTPrep_2DT_FW_MOCO_AVE_PSIR_STCNNT' 
    performUTValidation(UTCases, 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
    names{i}

    debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
    [data, gmap, res] = prepare_for_stcnnt_inference_LGE(debug_dir); 
    h = figure; imagescn(cat(4, data, res), [], [], [12], 3);
end

%% 3D

data_dir = '/data/raw_data/neuro/'

data_dir = '/data/raw_data/3D/'

data_dir = '/data/raw_data/noncardiac/neuro/'


data_dir = '/data/raw_data/noncardiac_20250805/neuro//'


UTCases = set_up_UT_cases_3D;

findAndMoveMeasDat(data_dir)
[names, num] = FindSubDirs(data_dir)

for i=1:num
    UTCases{1, 1} = data_dir;
    UTCases{1, 2} = names{i}

    h5name = fullfile(data_dir, names{i}, [names{i} '.h5']);
    dset = ismrmrd.Dataset(h5name);
    header = ismrmrd.xml.deserialize(dset.readxml());

    if header.encoding.encodingLimits.kspace_encoding_step_2.maximum > 2
        UTCases{1, 4} = 'GTPrep_3DT_Cartesian_AI_FIL.xml'
        UTCases{1, 5} = 'res_GTPrep_3DT_Cartesian_AI_FIL'

        UTCases{1, 4} = 'Generic_Cartesian_3D_Grappa_STCNNT.xml'
        UTCases{1, 5} = 'res_Generic_Cartesian_3D_Grappa_STCNNT'

    else
        UTCases{1, 4} = 'GTPrep_2DT_STCNNT_Spine.xml'
        UTCases{1, 5} = 'res_GTPrep_2DT_STCNNT_Slice'
    end
        

    UTCases{1, 3} = 'NX20'

    if strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0
            UTCases{1, 4} = 'default_measurement_dependencies.xml'
    end

    performUTValidation(UTCases(1,:), 0, 1, 'localhost', '9003', 1, 1, 0, 0, 0, [], [], 'wip_070_fire_IsmrmrdParameterMap_Siemens.xml', '/home/xueh/Debug/DebugOutput')
    names{i}

    % UTCases{1, 4} = 'Generic_RTCine.xml'
    % UTCases{1, 5} = 'res_Generic_RTCine' 

    debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
    [data, gmap, res] = prepare_for_stcnnt_inference_3D(debug_dir); 
    h = figure('Name', names{i}); imagescn(cat(4, data, res), [0 36*abs(mean(data(:)))], [size(data, 4), 2], [12], 3);
    saveas(h, fullfile(data_dir, [names{i} '.fig']));
    save(fullfile(data_dir, [names{i} '.mat']), 'data', 'gmap', 'res');
end

m_dir = 'model_res_73'

res_dir='/data/raw_data/noncardiac_20250805/neuro/meas_MID00604_FID124563_t1_mprage_cor_p2_iso_ipat4/res_Generic_Cartesian_3D_Grappa_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(res_dir, m_dir), 1);

res_dir='/data/raw_data/noncardiac_20250805/neuro/meas_MID00603_FID124562_t1_mprage_cor_p2_iso_ipat3/res_Generic_Cartesian_3D_Grappa_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(res_dir, m_dir), 1);

res_dir='/data/raw_data/noncardiac_20250805/neuro/meas_MID00602_FID124561_t1_mprage_cor_p2_iso/res_Generic_Cartesian_3D_Grappa_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(res_dir, m_dir), 1);

res_dir='/data/raw_data/noncardiac_20250805/neuro/meas_MID00608_FID124567_t2_swi_tra_p2_1_5mm/res_Generic_Cartesian_3D_Grappa_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(res_dir, m_dir), 1);

res_dir='/data/raw_data/noncardiac_20250805/neuro/meas_MID00598_FID124557_t2_tse_tra_512/res_GTPrep_2DT_STCNNT_Slice/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(res_dir, m_dir), 1);

res_dir='/data/raw_data/noncardiac_20250805/neuro/meas_MID00599_FID124558_t2_tse_tra_512_spat3/res_GTPrep_2DT_STCNNT_Slice/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(res_dir, m_dir), 1);

mat_dir = fullfile(data_dir, 'neuro_mat')
mkdir(mat_dir)

for i=3:num
    UTCases{1, 1} = data_dir;
    UTCases{1, 2} = names{i}

    h5name = fullfile(data_dir, names{i}, [names{i} '.h5']);
    dset = ismrmrd.Dataset(h5name);
    header = ismrmrd.xml.deserialize(dset.readxml());
    dset.close();
    
    if header.encoding.encodingLimits.kspace_encoding_step_2.maximum > 2
        UTCases{1, 4} = 'GTPrep_3DT_Cartesian_AI_FIL.xml'
        UTCases{1, 5} = 'res_GTPrep_3DT_Cartesian_AI_FIL'

        UTCases{1, 4} = 'Generic_Cartesian_3D_Grappa_STCNNT.xml'
        UTCases{1, 5} = 'res_Generic_Cartesian_3D_Grappa_STCNNT'

    else
        UTCases{1, 4} = 'GTPrep_2DT_STCNNT_Spine.xml'
        UTCases{1, 5} = 'res_GTPrep_2DT_STCNNT_Slice'
    end

    res_dir=fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput', m_dir);

    [input, res, gmap] = load_results_stcnnt_inference_perf(res_dir, 0);

    h = figure('Name', names{i});
    imagescn(cat(4, abs(input), abs(res), gmap), [], [1 3], [12], 3);
    save(fullfile(mat_dir, [names{i} '.mat']), 'input', 'res', 'gmap');
    saveas(h, fullfile(mat_dir, [names{i} '.fig']), 'fig');
end


%% LF

data_dir = '/data/raw_data/LF/FullExamExample/lge/'

findAndMoveMeasDat(data_dir)

UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;

[names, num] = FindSubDirs(data_dir)

for i=1:num
    UTCases{1, 1} = data_dir;
    UTCases{1, 2} = names{i}

    UTCases{1, 4} = 'Freemax_LGE_REP_STCNNT_offline.xml'
    UTCases{1, 5} = 'res_Freemax_LGE_REP_STCNNT_offline' 

    % UTCases{1, 4} = 'Freemax_LGE_REP_STCNNT_offline_half.xml'
    % UTCases{1, 5} = 'res_Freemax_LGE_REP_STCNNT_offline_half'
    % 
    UTCases{1, 4} = 'Freemax_LGE_REP_STCNNT.xml'
    UTCases{1, 5} = 'res_Freemax_LGE_REP_STCNNT' 
     
    if any([(strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0) (strfind(names{i}, 'noise_') > 0)])
            UTCases{1, 4} = 'default_measurement_dependencies_Noise_CoilSen_SCC.xml'
    end

    performUTValidation(UTCases, 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
    names{i}

    debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
    [data, gmap, res] = prepare_for_stcnnt_inference_LGE(debug_dir); 
    % if ~isempty(res)
    %     SLC = size(data, 5)
    %     if SLC==1
    %         h = figure('Name', names{i}); imagescn(cat(4, data, res), [], [SLC 4], [12], 3);
    %     else
    %         h = figure('Name', names{i}); imagescn(cat(4, data, res), [], [4 SLC], [4], 3);
    %     end
    %     saveas(h, fullfile(data_dir, [names{i} '_' UTCases{1, 5} '.fig']));
    % end
end

res_dir = '/data/raw_data/LGE_Denoising/Freemax_XL_NIH_120500_FID028851_G33_4CH_FB_de_tpat3_res256_Ave16_BW610_PHASres84/res_Freemax_LGE_REP_STCNNT_offline/DebugOutput//model_res'
[input, res, gmap] = load_results_stcnnt_inference_perf(res_dir, 1);

for i=1:num
    [a, header, acq_time, physio_time, endo_pt, epi_pt, user_int] =readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, 'res_Freemax_LGE_REP_STCNNT'), 111, 1, 0);
    % b = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, 'res_Freemax_LGE_REP_offline'), 111);
    % c = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, 'res_Freemax_LGE_REP_STCNNT_offline_half'), 111);

    disp([names{i} ' - ' num2str(size(a))])
    
    SLC = size(a, 3)
    for slc=1:SLC
        h = figure('name', names{i});
        %imagescn(cat(3, a(:,:,slc), c(:,:,slc), b(:,:,slc)), [header(slc).window_center-header(slc).window_width/2 header(slc).window_center+header(slc).window_width/2], [1 3], [12]);
        imagescn(a(:,:,slc), [header(slc).window_center-header(slc).window_width/2 header(slc).window_center+header(slc).window_width/2], [], [12]);
    
        saveas(h, fullfile(data_dir, [names{i} '_AI_slc' num2str(slc-1) '.fig']), 'fig');
    end
end

% -------------------------------------

data_dir = '/fastdata/lowfield//'

findAndMoveMeasDat(data_dir)

UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;

[names, num] = FindSubDirs(data_dir)

for i=1:num
    UTCases{1, 1} = data_dir;
    UTCases{1, 2} = names{i}

    UTCases{1, 4} = 'GT_RetroCine_AI.xml'
    UTCases{1, 5} = 'res_GT_RetroCine_AI' 
     
    if any([(strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0) (strfind(names{i}, 'noise_') > 0)])
            UTCases{1, 4} = 'default_measurement_dependencies_Noise_CoilSen_SCC.xml'
    end

    performUTValidation(UTCases(1,:), 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], 'wip_070_fire_IsmrmrdParameterMap_Siemens.xml', '/home/xueh/Debug/DebugOutput');
    names{i}

    % debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
    % [data, gmap, res] = prepare_for_stcnnt_inference_cine(debug_dir); 
    % if ~isempty(res)
    %     SLC = size(data, 4)
    %     if SLC< 3
    %         h = figure('Name', names{i}); imagescn(cat(4, data, res), [], [SLC 2], [12], 3);
    %     else
    %         h = figure('Name', names{i}); imagescn(cat(4, data, res), [], [4 ceil(SLC/2)], [8], 3);
    %     end
    %     saveas(h, fullfile(data_dir, [names{i} '_' UTCases{1, 5} '.fig']));
    % end
end

for i=1:num
    [a, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 1, 1);
    [b, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 2101, 1);

    size(a)

    h = figure('name', names{i});
    SLC = size(a, 4)
    if SLC < 2
        imagescn(cat(4, a, b), [], [SLC 2], [18], 3);
    else
        imagescn(cat(4, a, b), [], [2 SLC], [18], 3);
    end

    saveas(h, fullfile(data_dir, [names{i} '.fig']), 'fig');
end


% -------------------------------------

data_dir = '/fastdata/perfusion_h5/'

findAndMoveMeasDat(data_dir)

UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;

[names, num] = FindSubDirs(data_dir)

for i=1:num
    UTCases{1, 1} = data_dir;
    UTCases{1, 2} = names{i}

    UTCases{1, 4} = 'GT_QPerf_AI_STCNNT_noAIF.xml'
    UTCases{1, 5} = 'res_GT_QPerf_AI_STCNNT_noAIF' 
     
    if any([(strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0) (strfind(names{i}, 'noise_') > 0)])
            UTCases{1, 4} = 'default_measurement_dependencies_Noise_CoilSen_SCC.xml'
    end

    performUTValidation(UTCases, 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
    names{i}

    debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
    [data, gmap, res] = prepare_for_stcnnt_inference_perf(debug_dir); 
    if ~isempty(res)
        SLC = size(data, 4)
        if SLC< 3
            h = figure('Name', names{i}); imagescn(cat(4, data, res), [], [SLC 2], [12], 3);
        else
            h = figure('Name', names{i}); imagescn(cat(4, data, res), [], [4 ceil(SLC/2)], [8], 3);
        end
        saveas(h, fullfile(data_dir, [names{i} '_' UTCases{1, 5} '.fig']));
    end
end

keyframe = 15; 
level = 4; 
max_iter_num_pyramid_level = [100 100 100 100]; 
LocalCCR_sigmaArg = 2.0; 
BidirectionalReg = 1; 
dissimilarity_thres = 1e-6; 
DivergenceFreeReg = 0; 
KLT_regularization_hilbert_strength = [18 16 14 12 10]; 
KLT_regularization_minimal_hilbert_strength = 24; 
num_levels_moco_among_model = 5; 
regularization_scale_factor_moco_among_model = 1.0; 
dissimilarity_LocalCCR_sigmaArg_moco_among_model = 1.0; 
num_perf_for_PD = 3; 
DebugFolder = []; 
verbose = 0; 

A = abs(res(:,:,:,3));
[dx, dy, warped, dxInv, dyInv] = Matlab_gt_perfusion_model_moco(single(A), keyframe, level, max_iter_num_pyramid_level, LocalCCR_sigmaArg, BidirectionalReg, dissimilarity_thres, DivergenceFreeReg, KLT_regularization_hilbert_strength, KLT_regularization_minimal_hilbert_strength, num_levels_moco_among_model, regularization_scale_factor_moco_among_model, dissimilarity_LocalCCR_sigmaArg_moco_among_model, num_perf_for_PD, DebugFolder, verbose); 

B = Matlab_gt_apply_deformation_field_reg_2D_series(abs(input(:,:,:,3)), dx, dy);

figure; imagescn(cat(4, warped, B), [], [], [], 3);


% -------------------------------------

data_dir = '/data/raw_data/LF/non-cardiac/spine//'

findAndMoveMeasDat(data_dir)

UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;

[names, num] = FindSubDirs(data_dir)

for i=1:num
    UTCases{1, 1} = data_dir;
    UTCases{1, 2} = names{i}

    UTCases{1, 4} = 'GTPrep_2DT_STCNNT_Spine.xml'
    UTCases{1, 5} = 'res_GTPrep_2DT_STCNNT_Spine'
     
    if any([(strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0) (strfind(names{i}, 'noise_') > 0)])
            UTCases{1, 4} = 'default_measurement_dependencies_Noise_CoilSen_SCC.xml'
    end

    performUTValidation(UTCases, 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
    names{i}

    debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
    [data, gmap, res] = prepare_for_stcnnt_inference_spine(debug_dir); 
    if ~isempty(res)
        SLC = size(data, 3)
        h = figure('Name', names{i}); imagescn(abs(cat(4, data, res)), [], [2 size(data, 4)], [8], 3);
        saveas(h, fullfile(data_dir, [names{i} '_' UTCases{1, 5} '.fig']));
    end
end

res_dir = '/data/raw_data/LF/non-cardiac/spine/20201230_121929_meas_MID00103_FID16692_t1_tse_sag_t-spine/res_GTPrep_2DT_STCNNT_Spine/DebugOutput//model_res'

[input, res, gmap] = load_results_stcnnt_inference_perf(res_dir, 1);


% -------------------------------------

data_dir = '/data/raw_data/LF/non-cardiac/head//'

findAndMoveMeasDat(data_dir)

UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;

[names, num] = FindSubDirs(data_dir)

for i=1:num
    UTCases{1, 1} = data_dir;
    UTCases{1, 2} = names{i}

    UTCases{1, 4} = 'GTPrep_2DT_Cartesian_GFactor.xml'
    UTCases{1, 5} = 'res_GTPrep_2DT_Cartesian_GFactor'
     
    if any([(strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0) (strfind(names{i}, 'noise_') > 0)])
            UTCases{1, 4} = 'default_measurement_dependencies_Noise_CoilSen_SCC.xml'
    end

    performUTValidation(UTCases, 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
    names{i}

    debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
    [data, gmap, res] = prepare_for_stcnnt_inference_spine(debug_dir); 
    if ~isempty(res)
        SLC = size(data, 3)
        h = figure('Name', names{i}); imagescn(abs(cat(4, data, res)), [], [2 size(data, 4)], [8], 3);
        saveas(h, fullfile(data_dir, [names{i} '_' UTCases{1, 5} '.fig']));
    end
end

res_dir = '/data/raw_data/LF/non-cardiac/spine/20201230_121929_meas_MID00103_FID16692_t1_tse_sag_t-spine/res_GTPrep_2DT_STCNNT_Spine/DebugOutput//model_res'

res_dir = '/data/raw_data/neuro/Data_3DT_41837_2360858741_2360858746_574_20250114-174047/res_Generic_Cartesian_3D_Grappa_STCNNT/DebugOutput//model_res'
[input, res, gmap] = load_results_stcnnt_inference_perf(res_dir, 1);


%% OpenRecon
% ---------------------------------

host = 'localhost'
port = '9002'

host = 'GCRAZGDL4012.westus3.cloudapp.azure.com'
port = '9012'

data_dir_gt = '/data/raw_data/phantomdatfiles/snr2'

findAndMoveMeasDat(data_dir_gt)

[names, num] = FindSubDirs(data_dir_gt)

i = 1

UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;
UTCases{1, 1} = data_dir_gt;
UTCases{1, 2} = names{i}

UTCases{1, 4} = 'OR_Generic_Cartesian_Grappa_SNR.xml'
UTCases{1, 5} = 'res_OR_Generic_Cartesian_Grappa_SNR'

res_dir_OR = fullfile(data_dir_gt, names{i}, UTCases{1, 5})
debug_dir_OR = fullfile(res_dir, 'DebugOutput')

performUTValidation(UTCases(1,:), 0, 0, host, port, 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
names{i}

data_OR = readGTPlusExportImageSeries_Squeeze(res_dir_OR, 1);
figure; imagescn(data_OR, [], [], [6], 3);

std_map_OR = readGTPlusExportImageSeries_Squeeze(res_dir_OR, 3000);
figure; imagescn(std_map_OR)

% ---------------------------------------------

host = 'localhost'
port = '9002'

data_dir_gt = '//data/raw_data/OR/rtcine/'

findAndMoveMeasDat(data_dir_gt)

[names, num] = FindSubDirs(data_dir_gt)

i = 1

UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;
UTCases{1, 1} = data_dir_gt;
UTCases{1, 2} = names{i}

UTCases{1, 4} = 'OR_Generic_RTCine_STCNNT.xml'
UTCases{1, 5} = 'res_OR_Generic_RTCine_STCNNT_SPAT'

res_dir_OR = fullfile(data_dir_gt, names{i}, UTCases{1, 5})
debug_dir_OR = fullfile(res_dir_OR, 'DebugOutput')

performUTValidation(UTCases(1,:), 0, 0, host, port, 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
names{i}

i = 2

UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;
UTCases{1, 1} = data_dir_gt;
UTCases{1, 2} = names{i}

UTCases{1, 4} = 'Generic_RTCine_STCNNT.xml'
UTCases{1, 5} = 'res_Generic_RTCine_STCNNT'

res_dir = fullfile(data_dir_gt, names{i}, UTCases{1, 5})
debug_dir = fullfile(res_dir, 'DebugOutput')

performUTValidation(UTCases(1,:), 0, 0, host, port, 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
names{i}

data_OR = readGTPlusExportImageSeries_Squeeze(res_dir_OR, 1);
size(data_OR)
data = readGTPlusExportImageSeries_Squeeze(res_dir, 1);
size(data)
figure; imagescn(cat(4, data, data_OR), [0 4*mean(data(:))], [], [12], 3);

data_OR = readGTPlusExportImageSeries_Squeeze(res_dir_OR, 2101);
size(data_OR)
data = readGTPlusExportImageSeries_Squeeze(res_dir, 2101);
size(data)
figure; imagescn(cat(4, data, data_OR), [0 4*mean(data(:))], [], [12], 3);

% ---------------------------------------------

host = 'localhost'
port = '9002'

data_dir_gt = '//data/raw_data/OR/spine/'


findAndMoveMeasDat(data_dir_gt)

[names, num] = FindSubDirs(data_dir_gt)

i = 1

UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;
UTCases{1, 1} = data_dir_gt;
UTCases{1, 2} = names{i}

UTCases{1, 4} = 'OR_Generic_Cartesian_2D_Grappa_STCNNT.xml'
UTCases{1, 5} = 'res_OR_Generic_Cartesian_2D_Grappa_STCNNT'

res_dir_OR = fullfile(data_dir_gt, names{i}, UTCases{1, 5})
debug_dir_OR = fullfile(res_dir_OR, 'DebugOutput')

performUTValidation(UTCases(1,:), 0, 0, host, port, 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
names{i}

i = 2

UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;
UTCases{1, 1} = data_dir_gt;
UTCases{1, 2} = names{i}

UTCases{1, 4} = 'GTPrep_2DT_STCNNT_Spine.xml'
UTCases{1, 5} = 'res_GTPrep_2DT_STCNNT_Spine'

res_dir = fullfile(data_dir_gt, names{i}, UTCases{1, 5})
debug_dir = fullfile(res_dir, 'DebugOutput')

performUTValidation(UTCases(1,:), 0, 0, host, port, 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
names{i}

data_OR = readGTPlusExportImageSeries_Squeeze(res_dir_OR, 1);
size(data_OR)
data = readGTPlusExportImageSeries_Squeeze(res_dir, 1);
size(data)
figure; imagescn(cat(4, data, data_OR), [], [], [12], 3);

data_OR = readGTPlusExportImageSeries_Squeeze(res_dir_OR, 2101);
size(data_OR)
data = readGTPlusExportImageSeries_Squeeze(res_dir, 2101);
size(data)
figure; imagescn(cat(4, data, data_OR), [], [], [12], 3);

% ---------------------------------------------

%% Tyger
data_dir = '/data/raw_data/tyger/'

[names, num] = FindSubDirs(data_dir)

host = 'localhost'
port = '9002'

host = 'GCRSANDBOX113.redmond.corp.microsoft.com'
port = '9002'

host = 'GCRAZGDL1497.westus3.cloudapp.azure.com'
port = '9002'


% ---------------------------------

i = 3

UTCases = set_up_UT_cases_RTCine;
UTCases{1, 1} = data_dir;
UTCases{1, 2} = names{i}

UTCases{1, 4} = 'OR_Generic_Cartesian_Grappa_SNR.xml'
UTCases{1, 5} = 'res_Generic_Cartesian_Grappa_SNR'


UTCases{1, 4} = 'OR_SNR_Debug.xml'
UTCases{1, 5} = 'res_OR_SNR_Debug'

names{i}

performUTValidation(UTCases(1,:), 0, 0, host, port, 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
names{i}

gmap = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 20);
std_map = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 3000);

data = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 1);
res = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 300);
h = figure('Name', names{i}); imagescn(cat(4, data, res), [0 4*abs(mean(data(:)))], [size(data, 4), 2], [12], 3);
saveas(h, fullfile(data_dir, [names{i} '.fig']));

cd(fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput'))

cmaps = readGTPlusExportData('coil_map_encoding_0');
figure; imagescn(abs(squeeze(cmaps)))

kspace = readGTPlusExportData('data_encoding_0_1');
size(kspace)
figure; imagescn(abs(squeeze(data)), [], [8 4], [], 4)

figure; imagescn(data, [], [], [], 3);
figure; imagescn(std_map);

% ---------------------------------

i = 1

UTCases = set_up_UT_cases_RTCine;
UTCases{1, 1} = data_dir;
UTCases{1, 2} = names{i}

UTCases{1, 4} = 'OR_Generic_Cartesian_3D_Grappa_STCNNT.xml'
UTCases{1, 5} = 'res_Generic_Cartesian_3D_Grappa_STCNNT'

performUTValidation(UTCases(1,:), 0, 0, host, port, 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
names{i}

data = analyze75read(fullfile(data_dir, names{i}, UTCases{1, 5}, 'Generic_Cartesian_3D_Grappa_STCNNT.xml_SLC0_CON0_PHS0_REP0_SET0_AVE0_1_1.hdr'));
res = analyze75read(fullfile(data_dir, names{i}, UTCases{1, 5}, 'Generic_Cartesian_3D_Grappa_STCNNT.xml_SLC0_CON0_PHS0_REP0_SET0_AVE0_1_2101.hdr'));
h = figure('Name', names{i}); imagescn(cat(4, data, res), [0 4*abs(mean(data(:)))], [size(data, 4), 2], [12], 3);
saveas(h, fullfile(data_dir, [names{i} '.fig']));

% ---------------------------------

UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;

i=2
UTCases{1, 1} = data_dir;
UTCases{1, 2} = names{i}
UTCases{1, 4} = 'OR_Generic_RTCine_STCNNT.xml'
UTCases{1, 5} = 'res_Generic_RTCine_STCNNT' 
    
performUTValidation(UTCases(1,:), 0, 0, host, port, 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
names{i}

debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
[data, gmap, res] = prepare_for_stcnnt_inference_cine(debug_dir); 
h = figure('Name', names{i}); 
if size(data, 4) <= 2
    imagescn(cat(4, data, res), [0 2*abs(mean(data(:)))], [size(data, 4), 2], [12], 3);
else
    imagescn(cat(4, data, res), [0 2*abs(mean(data(:)))], [2 size(data, 4)], [12], 3);
end
saveas(h, fullfile(data_dir, [names{i} '.fig']));

% ---------------------------------

i = 4

UTCases = set_up_UT_cases_RTCine;
UTCases{1, 1} = data_dir;
UTCases{1, 2} = names{i}

UTCases{1, 4} = 'OR_Generic_Cartesian_Grappa_SNR.xml'
UTCases{1, 5} = 'res_Generic_Cartesian_Grappa_SNR'


UTCases{1, 4} = 'OR_SNR_Debug.xml'
UTCases{1, 5} = 'res_OR_SNR_Debug'

performUTValidation(UTCases(1,:), 0, 0, host, port, 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
names{i}

gmap = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 20);
std_map = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 3000);

data = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 1);
res = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, names{i}, UTCases{1, 5}), 300);
h = figure('Name', names{i}); imagescn(cat(4, data, res), [0 4*abs(mean(data(:)))], [size(data, 4), 2], [12], 3);
saveas(h, fullfile(data_dir, [names{i} '.fig']));

figure; imagescn(std_map)

%% test demo
data_dir = '/data/raw_data/SNR_PAPER/SNR_PAPER/R5/'

[names, num] = FindSubDirs(data_dir)

i = 1

UTCases{1, 1} = data_dir;
UTCases{1, 2} = names{i}
UTCases{1, 4} = 'Generic_Cartesian_OpenRecon_Demo.xml'
UTCases{1, 5} = 'res_Generic_Cartesian_OpenRecon_Demo'

performUTValidation(UTCases(1,:), 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
names{i}

debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput')

cd /data/raw_data/SNR_PAPER/SNR_PAPER/R5/RT_Cine_LIN_00000_226002092_226002101_4712_00000000-000000/res_Generic_Cartesian_OpenRecon_Demo/DebugOutput/

a = readGTPlusExportData('reports');
size(a)

figure; imagescn(abs(a))


a = readGTPlusExportData('rgb_2D');
size(a)

figure; imagescn(abs(squeeze(a)), [], [], [], 4)

mov = immovie(abs(a));
implay(mov);



a = readGTPlusExportData('rgb_3D');
size(a)

mov = immovie(permute(a, [1 2 4 3])/255.0);
implay(mov);

figure; imagescn(abs(squeeze(a)), [], [1 3], [], 3)


a = readGTPlusExportData('im_contours_landmarks');
size(a)

figure; imagescn(abs(squeeze(a)), [], [1 3], [], 3)


%% debug

cd(debug_dir)

a = analyze75read('PerfFlowMapping_completedJob_0_FpPreMap.hdr');
size(a)


b = analyze75read('PerfFlowMapping_completedJob_0_FlowMap.hdr');
size(b)

figure; imagescn(cat(3, a, b), [0 8], [], [8]); PerfColorMap;


x = readNPY('x.npy');
x2 = permute(x, [4 5 3 2 1]);
y = readNPY('y.npy');
y2 = permute(y, [4 5 3 2 1]);

figure; imagescn(x2(:,:,:,1,:), [], [], [], 3);
figure; imagescn(y2(:,:,:,1,:), [], [], [], 3);

d = readNPY('d.npy');
size(d)

d2 = permute(d, [2 3 1]);
size(d2)

figure; imagescn(d2, [], [], [1], 3);

t = readNPY('t.npy');
size(t)
t2 = permute(t, [2 3 1]);
size(t2)

figure; imagescn(t2, [], [], [1], 3);

%% test models

res_dir = 'model_res'
res_dir = 'model_res_4'
res_dir = 'model_res_5'
res_dir = 'model_res_6'
res_dir = 'model_res_7'
res_dir = 'model_res_8'
res_dir = 'model_res_9'

res_dir = 'model_res_10'

res_dir = 'model_res_11'

res_dir = 'model_res_12'

res_dir = 'model_res_13'

res_dir = 'model_res_14'

res_dir = 'model_res_15'

res_dir = 'model_res_16'

res_dir = 'model_res_18'

res_dir = 'model_res_19'

res_dir = 'model_res_20'

res_dir = 'model_res_29'

res_dir = 'model_res_30'

res_dir = 'model_res_31'

res_dir = 'model_res_32'

res_dir = 'model_res_33'

res_dir = 'model_res_35'

res_dir = 'model_res_36'

res_dir = 'model_res_37'

res_dir = 'model_res_38'

res_dir = 'model_res_39'

res_dir = 'model_res_40'

res_dir = 'model_res_41'


res_dir = 'model_res_42'

res_dir = 'model_res_43'

res_dir = 'model_res_44'

res_dir = 'model_res_45'

res_dir = 'model_res_46'

res_dir = 'model_res_47'

res_dir = 'model_res_48'

res_dir = 'model_res_49'

res_dir = 'model_res_50'

res_dir = 'model_res_51'

res_dir = 'model_res_52'

res_dir = 'model_res_53'

res_dir = 'model_res_54'

res_dir = 'model_res_55'

res_dir = 'model_res_56'

res_dir = 'model_res_57'

res_dir = 'model_res_58'

res_dir = 'model_res_59'

res_dir = 'model_res_60'

res_dir = 'model_res_70'

res_dir = 'model_res_76'

res_dir = 'model_res_77'

res_dir = 'model_res_79'

res_dir = 'model_res_80'

closeall
res_dir

data_dir = '/data/raw_data/perf_high_res/for_paper/Perfusion_AIF_Q_mapping_000000_441437177_441437186_821_00000000-000000/res_GT_QPerf_AI_STCNNT_OFFLINE/DebugOutput/for_model/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/noncardiac_20250805/neuro/meas_MID00599_FID124558_t2_tse_tra_512_spat3/res_GTPrep_2DT_STCNNT_Slice/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/noncardiac_20250805/knee/meas_MID00634_FID124593_t2_tse_sag_spat2_512/res_GTPrep_2DT_STCNNT_Spine/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir = '/data/raw_data/noncardiac_20250805/neuro/meas_MID00604_FID124563_t1_mprage_cor_p2_iso_ipat4/res_Generic_Cartesian_3D_Grappa_STCNNT/DebugOutput'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);


 
data_dir='/fastdata/lowfield/Freemax_XL_NIH_2024-11-13-122438_FID009204_G33_3ch_cine_256_R3ipat/res_GT_RetroCine_AI/DebugOutput'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/wb_psir/WB_LGE_MOCO_AVE_STCNNT_42537_749884520_749884527_280_20250428-121450/res_GTPrep_WB_LGE_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/perf_high_res/for_paper/Perfusion_AIF_Q_mapping_000000_441436899_441436908_363_00000000-000000/res_GT_QPerf_AI_STCNNT_OFFLINE/DebugOutput/for_model/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R5/RT_Cine_LIN_00000_226002173_226002182_4865_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R5/RT_Cine_LIN_00000_237143837_237143846_525_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R5/RT_Cine_LIN_00000_226002119_226002128_4761_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R6/RT_Cine_LIN_00000_237143644_237143653_178_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R6/RT_Cine_LIN_00000_226002119_226002128_4762_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R6/RT_Cine_LIN_00000_226002092_226002101_4713_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R6/RT_Cine_LIN_00000_226002119_226002128_4765_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R6/RT_Cine_LIN_00000_226002146_226002155_4815_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R6/RT_Cine_LIN_00000_226002173_226002182_4864_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R6/RT_Cine_LIN_00000_237143617_237143626_108_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R6/RT_Cine_LIN_00000_226002173_226002182_4864_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R6/RT_Cine_LIN_00000_226002146_226002155_4815_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R6/RT_Cine_LIN_00000_226002119_226002128_4765_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R6/RT_Cine_LIN_00000_226002092_226002101_4713_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);




data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R5/RT_Cine_LIN_00000_226002119_226002128_4764_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R5/RT_Cine_LIN_00000_226002146_226002155_4814_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R5/RT_Cine_LIN_00000_237143617_237143626_107_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R5/RT_Cine_LIN_00000_237143644_237143653_177_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R5/RT_Cine_LIN_00000_237143671_237143680_220_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R5/RT_Cine_LIN_00000_237143729_237143738_329_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/SNR_PAPER/SNR_PAPER/R5/RT_Cine_LIN_00000_226002092_226002101_4712_00000000-000000/res_Generic_RTCine_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/perf_high_res/for_paper/Perfusion_AIF_Q_mapping_000000_441437177_441437186_822_00000000-000000/res_GT_QPerf_AI_STCNNT_OFFLINE/DebugOutput/for_model/result'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);;

data_dir = '/data/raw_data/20250122_SPINE_DATFILES/meas_MID01532_FID37073__t2_tse_sag_p2_384/res_GTPrep_2DT_STCNNT_Spine/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir = '/data/raw_data/rtcine/RT_Cine_LIN_STCNNT_184270_442208992.2827355704.161462925_442208992.2827355704.161462925_1054_20250429-164003/res_Generic_RTCine_STCNNT/DebugOutput'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/rtcine/RT_Cine_LIN_STCNNT_184270_442208992.2827355704.161462925_442208992.2827355704.161462925_1053_20250429-163854/res_Generic_RTCine_STCNNT/DebugOutput/'




data_dir='/data/raw_data/neuro/Data_3DT_41837_2360858741_2360858746_573_20250114-173600/res_Generic_Cartesian_3D_Grappa_STCNNT/DebugOutput'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/neuro/meas_MID00073_FID33725_t1_mprage_sag_p3_p7mm/res_Generic_Cartesian_3D_Grappa_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/LGE_speedup/high_res/LGE_MOCO_AVE_STCNNT_41837_3717501238_3717501247_619_20241217-110806/res_GTPrep_WB_LGE_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/LF/non-cardiac/spine/20201230_120336_meas_MID00093_FID16682_t2_tse_sag_t-spine/res_GTPrep_2DT_STCNNT_Spine/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/LF/FullExamExample/cine/Freemax_XL_NIH_2024-11-13-122438_FID009204_G33_3ch_cine_256_R3ipat/res_Freemax_2DT_RetroCine_STCNNT_SCC/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/WB_LGE/good_cases/WB_LGE_MOCO_AVE_STCNNT_41837_2597033160_2597033169_332_20240212-164619/res_GTPrep_WB_LGE_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/LGE_speedup/regular/LGE_MOCO_AVE_OnTheFly_000000_479119192_479119199_920_00000000-000000/res_GTPrep_LGE_STCNNT/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/LF/perfusion_h5/h5/Freemax_XL_NIH_2024-11-13-120630_FID009184_G33_Adeno_stress_Perf_2RR_R3_192res_TI135/res_GT_QPerf_AI_STCNNT_noAIF/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/OR/spine/meas_MID01532_FID37073__t2_tse_sag_p2_384/res_GTPrep_2DT_STCNNT_Spine/DebugOutput//'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/perf/Perfusion_AIF_Q_mapping_202611_6.90.2.5001101287_6.90.2.5001101287_4166_20241213-095231/res_GT_QPerf_AI_STCNNT_OFFLINE/DebugOutput'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

data_dir='/data/raw_data/perf/Perfusion_AIF_Q_mapping_202611_6.90.2.5001105714_6.90.2.5001105714_199_20250321-092808/res_GT_QPerf_AI_STCNNT_OFFLINE/DebugOutput/'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);

%% test the kspace aliasing

cd /data/raw_data/retro_cine_3T/raw_selected/00047851/ch4/Retro_Lin_Cine_2DT_LAX_GLS_000000_591071194_591071203_877_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/

cd /data/raw_data/retro_cine_3T/raw_selected/00047995/sax/Retro_Lin_Cine_2DT_LAX_GLS_000000_613271213_613271222_297_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/

dst_dir = '/data/raw_data/retro_cine_3T/raw_selected'

[names, num] = FindSubDirs(dst_dir)

for i=385:num

    for k=1:3
        if k==1
            d_dir = fullfile(dst_dir, names{i}, 'ch2')
        end
        if k==2
            d_dir = fullfile(dst_dir, names{i}, 'ch4')
        end
        if k==3
            d_dir = fullfile(dst_dir, names{i}, 'sax')
        end

        [case_dirs, case_num] = FindSubDirs(d_dir)

        for n=1:case_num

            a_case_dir = fullfile(d_dir, case_dirs{n}, 'res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2', 'model');
            if ~exist(a_case_dir)
                continue
            end

            disp(['--> case ' num2str(i) ' - ' a_case_dir])

            cd(a_case_dir)

            input = complex(readNPY('input_real'), readNPY('input_imag'));
            size(input)
            
            data = complex(readNPY('data_real'), readNPY('data_imag'));
            size(data)
            
            ref = complex(readNPY('ref_real'), readNPY('ref_imag'));
            size(ref)
            
            unmix_c = complex(readNPY('unmix_c_real'), readNPY('unmix_c_imag'));
            size(unmix_c)
            
            coilMap = complex(readNPY('coil_maps_real'), readNPY('coil_maps_imag'));
            coilMap = double(squeeze(coilMap));
            size(coilMap)
            
            data = squeeze(data);
            size(data)
            
            ref = squeeze(ref);
            size(ref)
            
            RO = size(data, 1)
            E1 = size(data, 2)
            SLC = size(data, 5)
            PHS = size(data, 4)

            rec = load('for_model')
            
            fRO = readNPY('fRO.npy');
            norm(fRO(:))
            
            fE1 = readNPY('fE1.npy');
            norm(fE1(:))
            
            r_fRO = 1.0 / sqrt(sum(fRO.*fRO) / RO);
            r_fE1 = 1.0 / sqrt(sum(fE1.*fE1) / E1);            
            
            a = data(round(RO/2), :,1,1);
            ind = find(abs(a(:))>0);
            accelFactor_data = ind(2)-ind(1);
            
            kRO = 5;
            kE1 = 4;
            
            fitItself = 1;
            thres = [0.0005, 0.01 0.1 1.0 10.0];
            
            aug_dir = fullfile(a_case_dir, 'dealiasing');
            mkdir(aug_dir)

            for p=1:3
                kspace = data;
                accelFactor = accelFactor_data;
                note_str = 'spat';

                if p==2
                    for slc=1:SLC
                        for phs=1:2:PHS
                            kspace(:, ind(1:2:end), :, phs, slc) = 0;
                        end
                        for phs=2:2:PHS
                            kspace(:, ind(2:2:end), :, phs, slc) = 0;
                        end
                    end
                    
                    accelFactor = accelFactor_data* 2

                    note_str = 'tpat'
                end

                if p==3
                    for slc=1:SLC
                        kspace(:, ind(1:2:end), :, :, slc) = 0;
                    end
                    
                    accelFactor = accelFactor_data* 2

                    note_str = 'spat'
                end

                if p==4
                    for slc=1:SLC
                        kspace(:, ind(1:4:end), :, :, slc) = 0;
                    end
                    
                    accelFactor = accelFactor_data* 4

                    note_str = 'spat'
                end
            
                % compute the scaling ratio sqrt(N)/sqrt(W), W is the number of sampled
                % kspace points
                
                tt = kspace(:,:,1,1,1);
                W = numel(find(abs(100*tt(:))>0))
                
                snr_ratio = sqrt((RO*E1)/W)
                
                Im = [];
                gmaps = [];
                for k=1:numel(thres)
                    for slc=1:SLC
                        disp(['--> processing ' num2str(thres(k)) ' - ' num2str(slc)])
                        acsSrc = double(squeeze(ref(:,:,:,slc)));
                        acsDst = double(squeeze(ref(:,:,:,slc)));
                    
                        [ker, convKer] = Matlab_gt_grappa_2d_calibrate(acsSrc, acsDst, kRO, kE1, accelFactor, fitItself, thres(k));
                        
                        kIm = Matlab_gt_grappa_2d_compute_image_domain_kernel(convKer, RO, E1);
                        
                        [unmixing, gmap] = Matlab_gt_grappa_2d_compute_unmxing_coeff(kIm, coilMap(:,:,:,slc), accelFactor);
                               
                        a_im = snr_ratio * squeeze(sum(ifft2c(kspace(:,:,:,:,slc)) .* repmat(unmixing, [1 1 1 PHS]), 3));
                        gmaps(:,:,k,slc) = gmap;

                        filteredData = performKSpaceFilter2D(fft2c(a_im), fRO*r_fRO, fE1*r_fE1);

                        Im(:,:,:,k, slc) = ifft2c(filteredData);
                    end
                end
                size(Im)
                size(gmaps)

                writeNPY(real(Im), fullfile(aug_dir, ['im_dealiasing_R' num2str(accelFactor) '_' note_str '_real.npy']));
                writeNPY(imag(Im), fullfile(aug_dir, ['im_dealiasing_R' num2str(accelFactor) '_' note_str '_imag.npy']));
                writeNPY(gmaps, fullfile(aug_dir, ['gmaps_dealiasing_R' num2str(accelFactor) '_' note_str '.npy']));

                if p==1
                    im_clean = squeeze(Im(:,:,:,1,:));
                    gmap_clean = squeeze(gmaps(:,:,1,:));

                    writeNPY(real(im_clean), fullfile(aug_dir, ['im_dealiasing_R' num2str(accelFactor) '_clean_real.npy']));
                    writeNPY(imag(im_clean), fullfile(aug_dir, ['im_dealiasing_R' num2str(accelFactor) '_clean_imag.npy']));
                    writeNPY(gmap_clean, fullfile(aug_dir, ['gmaps_dealiasing_R' num2str(accelFactor) '_clean.npy']));
                end
            end
        end
    end
end


save Im Im gmaps

load('Im.mat')
size(Im)
size(gmaps)

figure; imagescn(abs(Im(:,:,:,:,4)), [], [], [], 3)

figure; imagescn(gmaps, [], [size(gmaps, 3), size(gmaps,4)],[12])

cd /data/raw_data/retro_cine_3T/raw_selected/00047848/ch2/Retro_Lin_Cine_2DT_LAX_GLS_000000_591071113_591071122_711_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/dealiasing/
cd /data/raw_data/retro_cine_3T/raw_selected/00047848/sax/Retro_Lin_Cine_2DT_LAX_GLS_000000_591071113_591071122_716_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/dealiasing

im = complex(readNPY('im_dealiasing_R4_tpat_real'), readNPY('im_dealiasing_R4_tpat_imag'));
size(im)
figure; imagescn(abs(im), [], [size(im, 4) size(im, 5)], [], 3)
gmaps = readNPY('gmaps_dealiasing_R4_tpat');
size(gmaps)
figure; imagescn(gmaps, [0 4], [size(im, 4) size(im, 5)])

im = complex(readNPY('im_dealiasing_R4_spat_real'), readNPY('im_dealiasing_R4_spat_imag'));
size(im)
figure; imagescn(abs(im), [], [size(im, 4) size(im, 5)], [], 3)
gmaps = readNPY('gmaps_dealiasing_R4_spat');
size(gmaps)
figure; imagescn(gmaps, [0 4], [size(im, 4) size(im, 5)])

im = complex(readNPY('im_dealiasing_R8_spat_real'), readNPY('im_dealiasing_R8_spat_imag'));
size(im)
figure; imagescn(abs(im), [], [size(im, 4) size(im, 5)], [], 3)
gmaps = readNPY('gmaps_dealiasing_R8_spat');
size(gmaps)
figure; imagescn(gmaps, [0 4], [size(im, 4) size(im, 5)])

im = complex(readNPY('im_dealiasing_R2_clean_real'), readNPY('im_dealiasing_R2_clean_imag'));
size(im)
figure; imagescn(abs(im), [], [], [], 3)
gmaps = readNPY('gmaps_dealiasing_R2_clean');
size(gmaps)
figure; imagescn(gmaps, [0 4])


%% clean the disk

[names, num] = FindSubDirs(dst_dir)

for i=1:num

    for k=1:3
        if k==1
            d_dir = fullfile(dst_dir, names{i}, 'ch2')
        end
        if k==2
            d_dir = fullfile(dst_dir, names{i}, 'ch4')
        end
        if k==3
            d_dir = fullfile(dst_dir, names{i}, 'sax')
        end

        [case_dirs, case_num] = FindSubDirs(d_dir)

        for n=1:case_num

            a_case_dir = fullfile(d_dir, case_dirs{n}, 'res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2', 'model', 'res');
            if ~exist(a_case_dir)
                continue
            end

            [slc_dirs, slc_num] = FindSubDirs(a_case_dir);
            for s=1:slc_num
                a_slc_dir = fullfile(a_case_dir, slc_dirs{s});
                [run_dirs, run_num] = FindSubDirs(a_slc_dir);
                for r=1:run_num
                    a_run_dir = fullfile(a_slc_dir, run_dirs{r})

                    try
                        if exist(readNPY(fullfile(a_run_dir, 'output_real.npy'))
                            res = complex(readNPY(fullfile(a_run_dir, 'output_real')), readNPY(fullfile(a_run_dir, 'output_imag')));
                            writeNPY(abs(res), fullfile(a_run_dir, 'output.npy'));
                            delete(fullfile(a_run_dir, 'output_real.npy'));
                            delete(fullfile(a_run_dir, 'output_imag.npy'));
                        end
                    catch
                        continue;
                    end
                end
            end
            
        end
    end
end


%% debug

cd /fastdata/mri/debug/

for k=0:9

    a = niftiread(['clean_im_' num2str(k) '.nii']);
    size(a)
    
    b = niftiread(['noisy_im_'  num2str(k) '.nii']);
    size(b)
    
    c = niftiread(['clean_im_degraded_'  num2str(k) '.nii']);
    size(c)
    
    d = niftiread(['gmap_'  num2str(k) '.nii']);
    size(d)
    
    
    h = figure; imagescn(cat(4, a, b, c, d), [], [1 4], [12], 3)
end

a = readNPY('/fastdata/mri/debug/data.npy');
b = readNPY('/fastdata/mri/debug/gmap.npy');
c = readNPY('/fastdata/mri/debug/noisy.npy');
d = readNPY('/fastdata/mri/debug/data_degraded.npy');

figure; imagescn(cat(3, a, b, c, d))

%% PhysioMRI

cd /data/raw_data/debug/Brain_Images/

% -----------------------------------------

a = load('brainIR_kspace3');

kspace = permute(a.kSpace3D, [2 3 1]);
size(kspace)

im = ifft3c(kspace);

mkdir data

writeNPY(single(real(im)), 'brainIR_real.npy');
writeNPY(single(imag(im)), 'brainIR_imag.npy');

gmap = ones(size(im));
writeNPY(single(gmap), 'brainIR_gmap.npy');

data_dir='/data/raw_data/debug/Brain_Images/data/'
res_dir = 'model_res_49'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);


% -----------------------------------------

a = load('brainT1_kspace');

kspace = permute(a.kSpace3D, [2 3 1]);
size(kspace)

im = ifft3c(kspace);
figure; imagescn(abs(im), [], [], [], 3)

data_dir='/data/raw_data/debug/Brain_Images/data/brainT1'
mkdir(data_dir)
cd(data_dir)

writeNPY(single(real(im)), 'input_real.npy');
writeNPY(single(imag(im)), 'input_imag.npy');

gmap = ones(size(im));
writeNPY(single(gmap), 'gmap.npy');

res_dir = 'model_res_49'
[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, res_dir), 1);


% -----------------------------------------

%% test data preparation

cd /data/raw_data/retro_cine_3T/raw_selected/00047855/ch4/Retro_Lin_Cine_2DT_LAX_GLS_000000_591071302_591071311_1097_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/dealiasing/

a = complex(readNPY('im_dealiasing_R4_spat_real'), readNPY('im_dealiasing_R4_spat_imag'));
size(a)

gmap = readNPY('gmaps_dealiasing_R4_spat');
size(gmap)

%% load in kspace data

cd /data/raw_data/wb_psir/WB_LGE_MOCO_AVE_STCNNT_41837_2597033160_2597033169_332_20240212-164619/res_GT_LGE/DebugOutput/

kspace = readGTPlusExportData('data_src_encoding_0');
size(kspace)

ref = readGTPlusExportData('ref_calib_encoding_0');
size(ref)

coil_map = readGTPlusExportData('coil_map_encoding_0');
size(coil_map)

im = readGTPlusExportData('unwrappedIm_encoding_0');
size(im)

%% compare data

cd /fastdata/data/realtime/RT_Cine_LIN_STCNNT_41837_118251232_118251241_2142_20250620-111517/res_Generic_RTCine_STCNNT/DebugOutput/
[data, gmap, res] = prepare_for_stcnnt_inference_cine('.');

cd /fastdata/data/realtime/RT_Cine_LIN_STCNNT_42537_72963109_72963116_69_20250620-072021/res_Generic_RTCine_STCNNT/DebugOutput/
[data2, gmap2, res2] = prepare_for_stcnnt_inference_cine('.');

figure; imagescn(gmap, [0 8]);colormap('jet')
figure; imagescn(gmap2, [0 8]);colormap('jet')

cd /fastdata/data/realtime
save two_cases data gmap res data2 gmap2 res2

load two_cases
figure; imagescn(gmap, [0 8]);colormap('jet')
figure; imagescn(cat(4, data, res, repmat(gmap, [1 1 size(data,3)])), [0 100], [1 3], [18], 3);


figure; imagescn(gmap2, [0 8]);colormap('jet')
figure; imagescn(cat(4, data2, res2, repmat(gmap2, [1 1 size(data2,3)])), [0 50], [1 3], [18], 3);


cd /data/raw_data/perf_high_res/for_paper/Perfusion_AIF_Q_mapping_000000_441437177_441437186_822_00000000-000000/res_GT_QPerf_AI_STCNNT_OFFLINE/DebugOutput/for_model/result/
a = complex(readNPY('input_real'), readNPY('input_imag'));
b = complex(readNPY('output_real'), readNPY('output_imag'));

figure; imagescn(abs(cat(4, a, b)), [], [2 3], [], 3);

cd /fastdata/data/wb_psir/WB_LGE_MOCO_AVE_OnTheFly_41837_1837798573_1837798582_338_20230703-134952/res_GTPrep_WB_LGE_STCNNT/DebugOutput/

psir = complex(analyze75read('PSIR_row0_moco_ave_slice_0_REAL'), analyze75read('PSIR_row0_moco_ave_slice_0_IMAG'));
psir_norm = complex(analyze75read('PSIR_Norm_row0_moco_ave_slice_0_REAL'), analyze75read('PSIR_Norm_row0_moco_ave_slice_0_IMAG'));
scc_map = analyze75read('sccMap_moco_ave_slice_0');

ir = complex(analyze75read('IRImage_moco_ave_slice_0_REAL'), analyze75read('IRImage_moco_ave_slice_0_IMAG'));
pd = complex(analyze75read('PDImage_moco_ave_slice_0_REAL'), analyze75read('PDImage_moco_ave_slice_0_IMAG'));

figure; imagescn(abs(cat(3, ir, pd)))

save test_data ir pd psir psir_norm scc_map

psir = ir .* conj(pd./abs(pd));
figure; imagescn(real(psir))

figure; imagescn(real(psir)./scc_map)


%% test

cd /home/xueh/mrprogs/imagingfm_BTCHW/projects/tests/data/snr_unit_data/

im = complex(readNPY('unwrappedIm_real'), readNPY('unwrappedIm_imag'));
size(im)

mask = readNPY('mask.npy');
size(mask)

im2 = readNPY('im_low_matrix');
size(im2)

mask2 = readNPY('mask_low_matrix.npy');
size(mask2)

%% check data

debug_dir = '/fastdata/mri/debug/'

closeall
for i=0:2
    clean_im = niftiread(fullfile(debug_dir, ['clean_im_' num2str(i) '.nii']));
    noisy_im = niftiread(fullfile(debug_dir, ['noisy_im_' num2str(i) '.nii']));
    gmap = niftiread(fullfile(debug_dir, ['gmap_' num2str(i) '.nii']));
    noise_sigmas = readNPY(fullfile(debug_dir, ['noise_sigmas_' num2str(i) '.npy']));

    size(clean_im)
    size(noisy_im)
    size(gmap)

    REP = size(clean_im, 4)
    for rep=1:REP
        clean_im(:,:,:,rep) = clean_im(:,:,:,rep) * noise_sigmas(rep);
        noisy_im(:,:,:,rep) = noisy_im(:,:,:,rep) * noise_sigmas(rep);
    end

    figure; imagescn(cat(4, clean_im, noisy_im, gmap), [], [3 REP], [12], 3);
end

%% check training

cd /home/xueh/mrprogs/resys-main/imaging-fm-projects/.run/output/train_samples

closeall
for b=0:10
    
    prefix = ['train_batch_' num2str(b)]

    noisy = complex(readNPY([prefix '_noisy_im_real']), readNPY([prefix '_noisy_im_imag']));
    size(noisy)
    
    clean = complex(readNPY([prefix '_clean_im_real']), readNPY([prefix '_clean_im_imag']));
    size(clean)
    
    gmap = readNPY([prefix '_gmap']);
    size(gmap)
    
    noise_sigma = readNPY([prefix '_noise_sigma']);
    noise_sigma

    h = figure; imagescn(cat(4, noisy, clean, gmap), [], [1 3], [], 3)

    try
        pred = complex(readNPY([prefix '_pred_im_real']), readNPY([prefix '_pred_im_imag']));
        size(pred)

        h = figure; imagescn(cat(4, noisy, clean, pred, gmap), [], [4 4], [], 3)
    catch
    end

    
end

cd /home/xueh/mrprogs/resys-main/imaging-fm-projects/.run/output

a = readNPY('image');
size(a)

b = readNPY('image_batch');
size(b)