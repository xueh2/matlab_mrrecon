
%% RetroCine

data_dir = '/home/gtuser/data1/retrocine/raw/'

findAndMoveMeasDat(data_dir)

UTCases = set_up_UT_cases_LGE;

[names, num] = FindSubDirs(data_dir)

names

for i=140:num

    disp([num2str(i) ' - ' names{i}])

    findAndMoveMeasDat(fullfile(data_dir, names{i}))
    [cases, K] = FindSubDirs(fullfile(data_dir, names{i}))

    
    for k=1:K
        UTCases{1, 1} = fullfile(data_dir, names{i});
        UTCases{1, 2} = cases{k};
    
        UTCases{1, 4} = 'GT_RetroCine_gmap_augmentation_OFFLINE.xml'
        UTCases{1, 5} = 'res_GT_RetroCine_gmap_augmentation_OFFLINE' 
    
        if strfind(cases{k}, 'ISMRMRD_Noise_dependency') > 0
            UTCases{1, 4} = 'default_measurement_dependencies_Noise_CoilSen_SCC.xml'
        end
    
        performUTValidation(UTCases, 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], [], '/home/gtuser/Debug/DebugOutput');
        cases{k}

    end

    delete('/tmp/gadgetron/*.h5')

    for k=1:K
        if strfind(cases{k}, 'ISMRMRD_Noise_dependency') > 0
            continue
        end

        debug_dir = fullfile(data_dir, names{i}, cases{k}, 'res_GT_RetroCine_gmap_augmentation_OFFLINE', 'DebugOutput');
    
        try
            [kspace, aliasedIm, acs_src, acs_dst, coil_map, recon_kspace, recon_im, ker, kIm, unmixing, gmap] = load_debug_output_for_recon_ai_retrocine(debug_dir, [2, 4], 0);
        
            py.numpy.save(fullfile(data_dir, names{i}, cases{k}, 'kspace.npy'), kspace);
            py.numpy.save(fullfile(data_dir, names{i}, cases{k}, 'aliasedIm.npy'), aliasedIm);
            py.numpy.save(fullfile(data_dir, names{i}, cases{k}, 'acs_src.npy'), acs_src);
            py.numpy.save(fullfile(data_dir, names{i}, cases{k}, 'acs_dst.npy'), acs_dst);
            py.numpy.save(fullfile(data_dir, names{i}, cases{k}, 'coil_map.npy'), coil_map);
            py.numpy.save(fullfile(data_dir, names{i}, cases{k}, 'recon_kspace.npy'), recon_kspace);
            py.numpy.save(fullfile(data_dir, names{i}, cases{k}, 'recon_im.npy'), recon_im);
            %py.numpy.save(fullfile(data_dir, names{i}, 'ker.npy'), ker);
            py.numpy.save(fullfile(data_dir, names{i}, cases{k}, 'kIm.npy'), kIm);
            py.numpy.save(fullfile(data_dir, names{i}, cases{k}, 'unmixing.npy'), unmixing);
            py.numpy.save(fullfile(data_dir, names{i}, cases{k}, 'gmap.npy'), gmap);
        
            command = ['azcopy copy ' fullfile(data_dir, names{i}, cases{k}) '  "https://mriimaging.blob.core.windows.net/data/mri/recon_ai_retro_cine/numpy/' names{i} '?${SAS}" --recursive  --include-pattern="*.npy"']
            dos(command ,'-echo')
        catch
        end

        try
            rmdir(fullfile(data_dir, names{i}, cases{k}, 'res_GT_RetroCine_gmap_augmentation_OFFLINE'), 's')
            delete(fullfile(data_dir, names{i}, cases{k}, '*.npy'));
        catch
        end
    end
end



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



cd /home/gtuser/data/raw/00055351/Perfusion_AIF_Q_mapping_000000_2491715859_2491715868_418_00000000-000000/res_GT_QPerf_AI_STCNNT_OFFLINE/DebugOutput

endo = analyze75read('Bullseye_endo_mask');
fmap = analyze75read('Bullseye_fmap');
rvi = analyze75read('Bullseye_rvi_mask');

size(endo)
size(fmap)
size(rvi)

[rvi_y, rvi_x] = find(rvi(:,:,1)>0)

[endo_y, endo_x] = find(endo(:,:,1)>0)

c_x = mean(endo_x)
c_y = mean(endo_y)

figure; imagescn(endo(:,:,1)); hold on; plot(c_x, c_y, '+');plot(rvi_x, rvi_y, 'o'); plot(12, 120, 'b*');; plot(0, 0, 'r*');

v1 = [rvi_x rvi_y] - [c_x c_y]
v2 = [cos(pi*240/180) sin(pi*240/180)]

v1 = v1 / norm(v1)

rot_deg = acosd(dot(v1,v2)/(norm(v1)*norm(v2)))
rot_deg = acos(dot(v1,v2)/(norm(v1)*norm(v2)))

% new rvi
res = [cos(rot_deg) -sin(rot_deg); sin(rot_deg) cos(rot_deg)] * [rvi_x-c_x; rvi_y-c_y]

rvi_x_n = res(1) + c_x
rvi_y_n = res(2) + c_y

figure; imagescn(endo(:,:,1)); hold on; plot(c_x, c_y, '+');plot(rvi_x, rvi_y, 'o'); plot(12, 120, 'b*');; plot(0, 0, 'r*'); % plot(rvi_x_n, rvi_y_n, 'y+', 'MarkerSize', 12);

% rotate mask
RO = size(endo, 1)
E1 = size(endo, 2)
new_endo_mask = zeros(RO, E1);
new_fmap = zeros(RO, E1);

[src_e1, src_ro] = meshgrid(1:E1, 1:RO);
I = endo(:,:,1);

I2 = fmap(:,:, 1);

n_c_x = RO/2
n_c_y = E1/2

n_c_x = c_x
n_c_y = c_y

for e1=1:E1
    for ro=1:RO

        % compute pixel in old images
        xq = cos(rot_deg) * (ro-1 - n_c_x) + sin(rot_deg) * (e1-1-n_c_y) + c_x;
        yq = -sin(rot_deg) * (ro-1 - n_c_x) + cos(rot_deg) * (e1-1-n_c_y) + c_y;

        if xq>0 & xq<RO & yq>0 & yq<E1
            new_endo_mask(ro, e1) = interp2(src_e1, src_ro, I, yq, xq, 'linear');
            new_fmap(ro, e1) = interp2(src_e1, src_ro, I2, yq, xq, 'linear');
        end

    end
end

new_endo_mask(find(isnan(new_endo_mask))) = 0;
new_fmap(find(isnan(new_endo_mask))) = 0;

figure; imagescn(new_endo_mask); hold on; plot(c_x, c_y, '+'); plot(rvi_x_n, rvi_y_n, 'o', 'MarkerSize', 12, 'LineWidth', 3);


figure; imagescn(new_fmap, [0 8]); PerfColorMap; hold on; plot(c_x, c_y, '+'); plot(rvi_x_n, rvi_y_n, 'go', 'MarkerSize', 12, 'LineWidth', 3);

%% get neuro and spine figures

cd /data/raw_data/noncardiac/noncardiac/case5_res/Data_3DT_66016_056941164_056941169_33_20251023-112717/model_res_81/

[input1, res1, gmap1] = load_results_stcnnt_inference_perf('/home/gtuser/data/noncardiac/case5_res/Data_3DT_66016_056941164_056941169_33_20251023-112717/', 0);
[input2, res2, gmap2] = load_results_stcnnt_inference_perf('/home/gtuser/data/noncardiac/case5_res/Data_3DT_66016_056941164_056941169_34_20251023-113309/', 0);
[input3, res3, gmap3] = load_results_stcnnt_inference_perf('/home/gtuser/data/noncardiac/case5_res/Data_3DT_66016_056941164_056941169_35_20251023-113656/', 0);

size(res1)
size(res3)

a1 = interpimages3D(input1, 'linear', size(res3));
a2 = interpimages3D(input2, 'linear', size(res3));

out1 = interpimages3D(res1, 'linear', size(res3));
out2 = interpimages3D(res2, 'linear', size(res3));

h = figure; imagescn(cat(4, a1, a2, input3, out1, out2, res3), [], [2 3], [16], 3)
saveas(h, '/data/raw_data/noncardiac/noncardiac/case5_brain_T1.fig', 'fig')

[input1, res1, gmap1] = load_results_stcnnt_inference_perf('/data/raw_data/noncardiac/noncardiac/case5_res/Data_2D_66016_056941164_056941169_30_20251023-112121/model_res_81/', 0);
[input2, res2, gmap2] = load_results_stcnnt_inference_perf('/data/raw_data/noncardiac/noncardiac/case5_res/Data_2D_66016_056941164_056941169_31_20251023-112322/model_res_81/', 0);
[input3, res3, gmap3] = load_results_stcnnt_inference_perf('/data/raw_data/noncardiac/noncardiac/case5_res/Data_2D_66016_056941164_056941169_32_20251023-112431/model_res_81/', 0);

[input1, res1, gmap1] = load_results_stcnnt_inference_perf('/data/raw_data/noncardiac/noncardiac/case5_res/Spine_2D_66016_056941164_056941169_51_20251023-115001/model_res_81/', 0);
[input2, res2, gmap2] = load_results_stcnnt_inference_perf('/data/raw_data/noncardiac/noncardiac/case5_res/Spine_2D_66016_056941164_056941169_52_20251023-115206/model_res_81/', 0);
[input3, res3, gmap3] = load_results_stcnnt_inference_perf('/data/raw_data/noncardiac/noncardiac/case5_res/Spine_2D_66016_056941164_056941169_53_20251023-115340/model_res_81/', 0);





