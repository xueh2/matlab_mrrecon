clear all
closeall

addpath('/data/raw_data/SNR_PAPER/')

%% prepare 3T for Cine contouring

data_dir = '/data/raw_data/retro_cine_3T/raw/'

dst_dir = '/data/raw_data/retro_cine_3T/raw_selected/'
mkdir(dst_dir)

% [raw_dirs, num_raw] = FindSubDirs(dst_dir);
% for k=1:num_raw-1
%     [sub_dirs, num_subs] = FindSubDirs(fullfile(dst_dir, raw_dirs{k}));
%     num_subs
%     if num_subs == 1
%         rmdir(fullfile(dst_dir, raw_dirs{k}), 's');
%     end
% end

dst_fig_dir = fullfile(dst_dir, 'figures')
mkdir(dst_fig_dir)

snr_levels = [0.05, 0.1, 0.2, 0.5, 0.75, 1.0, 1.5, 2.0, 4.0, 8.0]
num_levels = numel(snr_levels)

[names, num] = FindSubDirs(data_dir)

for n=1:num

    case_dir = fullfile(data_dir, names{n})
    [h5names, num_h5] = findFILE(case_dir, 'ISMRMRD_Noise_dependency*.h5');

    case_dst_dir = fullfile(dst_dir, names{n}, 'depen');
    mkdir(case_dst_dir)

    for t=1:num_h5
        [p, h5_name, ext] = fileparts(h5names{t});
        movefile(h5names{t}, case_dst_dir);
    end
    findAndMoveMeasDat(case_dst_dir);
end

for n=1:num

    case_dir = fullfile(data_dir, names{n})
    [h5names, num_h5] = findFILE(case_dir, 'Retro_Lin_Cine*.h5');

    if(num_h5>3)

        ch2 = [];
        ch4 = [];
        sax = [];

        for h=1:num_h5
            try
                dset = ismrmrd.Dataset(h5names{h});
                hdr = ismrmrd.xml.deserialize(dset.readxml);
                dset.close();
                pro_name = hdr.measurementInformation.protocolName;
                pro_name = lower(pro_name);
                if ~isempty(strfind(pro_name, '2ch')) | ~isempty(strfind(pro_name, 'ch2')) | ~isempty(strfind(pro_name, '2 ch')) | ~isempty(strfind(pro_name, 'ch 2')) | ~isempty(strfind(pro_name, 'vla'))
                    ch2 = [ch2; {h5names{h}}];
                end
    
                if ~isempty(strfind(pro_name, '4ch')) | ~isempty(strfind(pro_name, 'ch4')) | ~isempty(strfind(pro_name, '4 ch')) | ~isempty(strfind(pro_name, 'ch 4')) | ~isempty(strfind(pro_name, 'hla'))
                    ch4 = [ch4; {h5names{h}}];
                end
    
                if ~isempty(strfind(pro_name, 'sax')) | ~isempty(strfind(pro_name, 'sa'))
                    sax = [sax; {h5names{h}}];
            end
            catch
                continue
            end
        end

        disp([names{n} ' - ch2 ' num2str(numel(ch2)) ' - ch4 ' num2str(numel(ch4)) ' - sax ' num2str(numel(sax))])

        if numel(ch2)>0 & numel(ch4)>0 & numel(sax)>0
            case_dst_dir = fullfile(dst_dir, names{n});
            mkdir(case_dst_dir)

            d_dir = fullfile(case_dst_dir, 'ch2');
            mkdir(d_dir)
            for k=1:numel(ch2)
                movefile(ch2{k}, d_dir);
            end

            d_dir = fullfile(case_dst_dir, 'ch4');
            mkdir(d_dir)
            for k=1:numel(ch4)
                movefile(ch4{k}, d_dir);
            end

            d_dir = fullfile(case_dst_dir, 'sax');
            mkdir(d_dir)
            for k=1:numel(sax)
                movefile(sax{k}, d_dir);
            end
    
        end

    end
end

UTCases = set_up_UT_cases_Cine
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

        [h5names, num_h5] = findFILE(d_dir, 'Retro_Lin_Cine*.h5');
    
        for n=1:num_h5
            UTCases{1, 1} = d_dir;
            [p, h5name, ext] = fileparts(h5names{n});
            mkdir(fullfile(d_dir, h5name))
            movefile(h5names{n}, fullfile(d_dir, h5name));
        end
    end
end

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

        [h5_dirs, num_dirs] = FindSubDirs(d_dir);

        for n=1:num_dirs
            debug_dir = fullfile(d_dir, h5_dirs{n}, 'res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE');
            if exist(debug_dir)
                disp(['remove - ' debug_dir])
                rmdir(debug_dir, 's')            
            end
        end
    end
end

% run depen scans
for i=1:num

    d_dir = fullfile(dst_dir, names{i}, 'depen');
    [h5_dirs, num_dirs] = FindSubDirs(d_dir);

    for n=1:num_dirs

        UTCases{1, 1} = d_dir;
        UTCases{1, 2} = h5_dirs{n}

        UTCases{1, 4} = 'default_measurement_dependencies_Noise_CoilSen_SCC.xml'
        UTCases{1, 5} = 'res' 

        performUTValidation(UTCases, 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
        names{i}
    end
end

fig_picked_dir = fullfile('/data/raw_data/retro_cine_3T/figure_selected/');
mkdir(fig_picked_dir)
[fignames, fig_num] = findFILE(fig_picked_dir, '*sax*.fig')

xml_name = 'GTPrep_2DT_RetroCine_SCC_for_AI_OFFLINE.xml'
res_name = 'res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE' 
im_series_num = 0;
im_scc_series_num = 1;

xml_name = 'GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE.xml'
res_name = 'res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2' 
im_series_num = 1;
im_scc_series_num = 2;
gmaps_series_num = 300;

for i=1:fig_num

    [p, fname, ext] = fileparts(fignames{i});

    ind = find(fname=='_');

    case_name = fname(1:ind(1)-1)
  
    d_dir = fullfile(dst_dir, case_name, 'depen');
    [h5_dirs, num_depen_dirs] = FindSubDirs(d_dir);

    % for n=1:num_depen_dirs
    % 
    %     UTCases{1, 1} = d_dir;
    %     UTCases{1, 2} = h5_dirs{n}
    % 
    %     UTCases{1, 4} = 'default_measurement_dependencies_Noise_CoilSen_SCC.xml'
    %     UTCases{1, 5} = 'res' 
    % 
    %     performUTValidation(UTCases(1, :), 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
    %     names{i}
    % end

    for k=1:3
        if k==1
            d_dir = fullfile(dst_dir, case_name, 'ch2')
        end
        if k==2
            d_dir = fullfile(dst_dir,case_name, 'ch4')
        end
        if k==3
            d_dir = fullfile(dst_dir,case_name, 'sax')
        end

        [h5_dirs, num_dirs] = FindSubDirs(d_dir);
    
        for n=1:num_dirs
            UTCases{1, 1} = d_dir;
            UTCases{1, 2} = h5_dirs{n};
        
            UTCases{1, 4} = xml_name;
            UTCases{1, 5} = res_name;
        
            res_dir = fullfile(d_dir, h5_dirs{n}, UTCases{1, 5})
            model_dir = fullfile(res_dir, 'model');
            if exist(fullfile(model_dir, 'fE1.npy'))
                delete(fullfile(res_dir, '*.img'));
                delete(fullfile(res_dir, '*.attrib'));
                delete(fullfile(res_dir, '*.hdr'));
                delete(fullfile(res_dir, '*.h5'));
                continue
            end

            performUTValidation(UTCases(1, :), 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput');
            case_name            

       
            h5_fname = fullfile(d_dir, h5_dirs{n}, [h5_dirs{n} '.h5']);

            dset = ismrmrd.Dataset(h5_fname);
            hdr = ismrmrd.xml.deserialize(dset.readxml);
            dset.close();

            SLC = hdr.encoding.encodingLimits.slice.maximum+1;
            
            debug_dir = fullfile(res_dir, 'DebugOutput');

            try
                %[data, header, acq_time, physio_time, attribs, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(res_dir, im_series_num, 1, 0); size(data)
                %[data_scc, header_scc, acq_time, physio_time, attribs, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(res_dir, im_scc_series_num, 1, 0); size(data_scc)
                %[data, header, attribs, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_h5(folderName, seriesNum, withTime, numAsRep);

                [data, header, attribs, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_h5(res_dir, im_series_num, 1, 0); size(data)
                [gmaps, header_gmaps, attribs_gmaps, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_h5(res_dir, gmaps_series_num, 1, 0); size(gmaps)

                data = squeeze(data);
                gmaps = squeeze(gmaps);

                if SLC > 1
                    data = permute(data, [2 1 4 3]); size(data)
                    gmaps = permute(gmaps, [2 1 4 3]); size(gmaps)
                else
                    data = permute(data, [2 1 3]); size(data)
                    gmaps = permute(gmaps, [2 1 3]); size(gmaps)
                end

                %gmaps = readGTPlusExportImageSeries_Squeeze(res_dir, 300); size(gmaps)

                %gmaps_aug = analyze75read(fullfile(debug_dir, 'gfactor_augmented_encoding_0_1'));
                %size(gmaps_aug)
                gmaps_aug = gmaps;

                coil_maps = readGTPlusExportData(fullfile(debug_dir, 'coil_map_encoding_0')); size(coil_maps)

                ref_src = readGTPlusExportData(fullfile(debug_dir, 'ref_encoding_0_1')); size(ref_src)
                ref_dst = readGTPlusExportData(fullfile(debug_dir, 'ref_calib_dst_encoding_0_1')); size(ref_dst)
                norm(ref_src(:)-ref_dst(:))

                ref_for_coil_map = readGTPlusExportData(fullfile(debug_dir, 'ref_coil_map_encoding_0_1')); size(ref_for_coil_map)

                data_src = readGTPlusExportData(fullfile(debug_dir, 'data_src_encoding_0')); 
                size(data_src)

                RO = size(data_src, 1)
                E1 = size(data_src, 2)
                CHA = size(data_src, 4)

                unmix_c = [];
                for slc=1:SLC
                    for R=2:8
                        a_unmix_c = readGTPlusExportData(fullfile(debug_dir, ['GFactor_aug_unmixC_n0_s0_slc' num2str(slc-1) '_encoding_0_R' num2str(R)]));
                        if isempty(unmix_c)
                            unmix_c = zeros(RO, E1, CHA, SLC, 7);
                        end
                        unmix_c(:,:,:,slc, R-1) = a_unmix_c;
                    end
                end
                size(unmix_c)

                fRO = readGTPlusExportData(fullfile(debug_dir, 'filterRO_encoding_0_Image'));
                fE1 = readGTPlusExportData(fullfile(debug_dir, 'filterE1_encoding_0_Image'));

                size(fRO)
                size(fE1)                

            catch
                continue
            end

            RO = size(data, 1);
            E1 = size(data, 2);

            model_dir = fullfile(res_dir, 'model');
            mkdir(model_dir)

            if num_depen_dirs>0
                for slc=1:SLC
                    scc_map_file = fullfile(debug_dir, ['scc_map_slc_' num2str(slc-1) '.hdr']);
                    if ~exist(scc_map_file)
                        continue
                    end
                    a_scc = analyze75read(scc_map_file);
                    disp(['write out scc map for slice ' num2str(slc-1) ' ... '])
                    writeNPY(permute(a_scc, [2 1]), fullfile(model_dir, ['scc_maps_' num2str(slc-1) '.npy']));
                end
            end

            data = data / 10.0;
            gmaps = gmaps / 100.0;

           
            save(fullfile(model_dir, 'for_model.mat'), 'header', 'acq_time', 'physio_time', 'attribs');

            % data = permute(data, [2 1 4 3]);
            % size(data)
            % 
            % gmaps = permute(gmaps, [2 1 4 3]);
            % size(gmaps)
            % 
            % gmaps_aug = permute(gmaps_aug, [2 1 4 3]);
            % size(gmaps_aug)

            data = squeeze(data);
            gmaps = squeeze(gmaps);

            writeNPY(real(data), fullfile(model_dir, 'input_real.npy'));
            writeNPY(imag(data), fullfile(model_dir, 'input_imag.npy'));
            %writeNPY(gmaps(:,:,:,1), fullfile(model_dir, 'gmap.npy'));
            %writeNPY(squeeze(gmaps_aug), fullfile(model_dir, 'gmaps_aug.npy'));            

            writeNPY(squeeze(gmaps), fullfile(model_dir, 'gmaps_aug.npy'));

            writeNPY(real(coil_maps), fullfile(model_dir, 'coil_maps_real.npy'));
            writeNPY(imag(coil_maps), fullfile(model_dir, 'coil_maps_imag.npy'));

            writeNPY(real(ref_src), fullfile(model_dir, 'ref_real.npy'));
            writeNPY(imag(ref_src), fullfile(model_dir, 'ref_imag.npy'));

            writeNPY(real(data_src), fullfile(model_dir, 'data_real.npy'));
            writeNPY(imag(data_src), fullfile(model_dir, 'data_imag.npy'));

            writeNPY(real(unmix_c), fullfile(model_dir, 'unmix_c_real.npy'));
            writeNPY(imag(unmix_c), fullfile(model_dir, 'unmix_c_imag.npy'));

            writeNPY(abs(fRO), fullfile(model_dir, 'fRO.npy'));
            writeNPY(abs(fE1), fullfile(model_dir, 'fE1.npy'));

            rmdir(debug_dir, 's');
            delete(fullfile(res_dir, '*.img'));
            delete(fullfile(res_dir, '*.attrib'));
            delete(fullfile(res_dir, '*.hdr'));
            delete(fullfile(res_dir, '*.h5'));

            disp([fullfile(dst_dir, case_name)]);
        end
    end
end

% make figures
for i=1:num  
    disp([num2str(i) ' - ' names{i}])
    for k=1:3
        if k==1
            d_dir_s = 'ch2'
        end
        if k==2
            d_dir_s = 'ch4'
        end
        if k==3
            d_dir_s = 'sax'            
        end

        d_dir = fullfile(dst_dir, names{i}, d_dir_s)
        [h5_dirs, num_dirs] = FindSubDirs(d_dir);
    
        for n=1:num_dirs
            UTCases{1, 1} = d_dir;
            UTCases{1, 2} = h5_dirs{n};
        
            UTCases{1, 4} = xml_name;
            UTCases{1, 5} = res_name;
        
            res_dir = fullfile(d_dir, h5_dirs{n}, UTCases{1, 5});
            model_dir = fullfile(res_dir, 'model');

            if exist(fullfile(dst_fig_dir, [names{i} '_' d_dir_s '_' num2str(n) '.fig']))
                continue
            end

            try

                data = complex(readNPY(fullfile(model_dir, 'input_real')), readNPY(fullfile(model_dir, 'input_imag')));
                size(data)
    
                data = squeeze(data);
    
                mag_scc = abs(data);           
                size(mag_scc)
    
                SLC = size(data, 4);
    
                h = figure('Name', [d_dir_s '_' num2str(n) '_' names{i}]);
                if SLC == 1
                    imagescn(mag_scc, [0 4*mean(mag_scc(:))], [1 SLC], [12], 3);
                else if (SLC <= 3)
                    imagescn(mag_scc, [0 4*mean(mag_scc(:))], [1 SLC], [12], 3);
                else 
                    if SLC <=6
                        imagescn(mag_scc, [0 4*mean(mag_scc(:))], [2 ceil(SLC/2)], [12], 3);
                    else
                        imagescn(mag_scc, [0 4*mean(mag_scc(:))], [3 ceil(SLC/3)], [12], 3);
                    end
                end
                end
                saveas(h, fullfile(dst_fig_dir, [names{i} '_' d_dir_s '_' num2str(n) '.jpg']), 'jpg')
                saveas(h, fullfile(dst_fig_dir, [names{i} '_' d_dir_s '_' num2str(n)]), 'fig')
                close(h)
            catch
                continue
            end
        end
    end
end

% pick the cases

for i=1:num
  
    disp([num2str(i) ' - ' names{i}])


    ch2_fig = fullfile(dst_fig_dir, [names{i} '_ch2_1.fig']);    
    ch4_fig = fullfile(dst_fig_dir, [names{i} '_ch4_1.fig']);
    sax_fig = fullfile(dst_fig_dir, [names{i} '_sax_1.fig']);
    
    if exist(ch2_fig) & exist(ch4_fig) & exist(sax_fig)

        found_flag = 0;
        for fg=1:fig_num
            [fpath, fname, ext] = fileparts(fignames{fg});
            if strcmp(fname, [names{i} '_sax_1'])==1
                found_flag = 1;
                break;
            end
        end
    
        if found_flag
            continue
        end

        h1 = open(ch2_fig); 
        figure(h1);MR_toolbox;
        
        h2= open(ch4_fig); 
        figure(h2);MR_toolbox;
        
        h3= open(sax_fig); 
        figure(h3);MR_toolbox;

        reply = input('Do you want keep this case? Y/N [Y]:','s');
        if isempty(reply)
          reply = 'Y';
        end

        if reply == 'Y'
            movefile(ch2_fig, fig_picked_dir);
            movefile(ch4_fig, fig_picked_dir);
            movefile(sax_fig, fig_picked_dir);
        end
        closeall
    end
end

% assemble the training
[fignames, fig_num] = findFILE(fig_picked_dir, '*sax*.fig')

cases_list = [];
for n=1:fig_num

    [p, fname, ext] = fileparts(fignames{n});

    ind = find(fname=='_');

    case_name = fname(1:ind(1)-1)


    d_dir = fullfile(dst_dir, case_name, 'ch2');
    [h5_dirs, num_dirs] = FindSubDirs(d_dir);
    ch2_dir = fullfile(d_dir, h5_dirs{1})

    d_dir = fullfile(dst_dir, case_name, 'ch4');
    [h5_dirs, num_dirs] = FindSubDirs(d_dir);
    ch4_dir = fullfile(d_dir, h5_dirs{1})

    d_dir = fullfile(dst_dir, case_name, 'sax');
    [h5_dirs, num_dirs] = FindSubDirs(d_dir);
    sax_dir = fullfile(d_dir, h5_dirs{1})

    ch2_model_dir = fullfile(ch2_dir, res_name, 'model');
    ch4_model_dir = fullfile(ch4_dir, res_name, 'model');
    sax_model_dir = fullfile(sax_dir, res_name, 'model');

    if exist(fullfile(ch2_model_dir, 'gmaps_aug.npy'))
        cases_list = [cases_list; {ch2_model_dir}];
    end
    if exist(fullfile(ch4_model_dir, 'gmaps_aug.npy'))
        cases_list = [cases_list; {ch4_model_dir}];
    end
    if exist(fullfile(sax_model_dir, 'gmaps_aug.npy'))
        cases_list = [cases_list; {sax_model_dir}];
    end
end

size(cases_list)

R = [2, 3, 4, 5];

ignore_check_existed = 1;

for wg=0:3    
    which_gmap = wg
    run_file_name = fullfile(fig_picked_dir, ['create_snr_data_' res_name '_which_gmap_' num2str(which_gmap) '.sh'])
    create_a_test_case_snr(run_file_name, cases_list, which_gmap, ignore_check_existed, ['res_R_' num2str(R(which_gmap+1))])
end

% for i=1:size(cases_list, 1)
%     a_case = cases_list{i, 1}
%     rmdir(fullfile(a_case, 'res_R_2'), 's');
%     rmdir(fullfile(a_case, 'res_R_3'), 's');
%     rmdir(fullfile(a_case, 'res_R_4'), 's');
%     rmdir(fullfile(a_case, 'res_R_5'), 's');
% end

% get the case to run models
cases_list_ai = [];

res_dir_list = {'res'}

res_dir_list = {'res_R_2', 'res_R_3', 'res_R_4', 'res_R_5'}

for n=1:fig_num

    [p, fname, ext] = fileparts(fignames{n});

    ind = find(fname=='_');

    case_name = fname(1:ind(1)-1)


    d_dir = fullfile(dst_dir, case_name, 'ch2');
    [h5_dirs, num_dirs] = FindSubDirs(d_dir);
    ch2_dir = fullfile(d_dir, h5_dirs{1})

    d_dir = fullfile(dst_dir, case_name, 'ch4');
    [h5_dirs, num_dirs] = FindSubDirs(d_dir);
    ch4_dir = fullfile(d_dir, h5_dirs{1})

    d_dir = fullfile(dst_dir, case_name, 'sax');
    [h5_dirs, num_dirs] = FindSubDirs(d_dir);
    sax_dir = fullfile(d_dir, h5_dirs{1})

    for p=1:numel(res_dir_list)
        a_res = res_dir_list{p};
        ch2_model_dir = fullfile(ch2_dir, res_name, 'model', a_res);
        ch4_model_dir = fullfile(ch4_dir, res_name, 'model', a_res);
        sax_model_dir = fullfile(sax_dir, res_name, 'model', a_res);
    
        [f, nn] = FindSubDirs(ch2_model_dir);
        for k=1:nn
            cases_list_ai = [cases_list_ai; {fullfile(ch2_model_dir, f{k})}];
        end
        [f, nn] = FindSubDirs(ch4_model_dir);
        for k=1:nn
            cases_list_ai = [cases_list_ai; {fullfile(ch4_model_dir, f{k})}];
        end
        [f, nn] = FindSubDirs(sax_model_dir);
        for k=1:nn
            cases_list_ai = [cases_list_ai; {fullfile(sax_model_dir, f{k})}];
        end
    end
end

size(cases_list_ai)

% for kk=1:numel(cases_list_ai)
%     a_case = cases_list_ai{kk};
%     gmaps = readNPY(fullfile(a_case, "../../gmaps_aug.npy"));
%     disp([a_case ' - ' num2str(size(gmaps))])
%     writeNPY(squeeze(gmaps(:,:,1,:)), fullfile(a_case, "../../gmap.npy"));
% end

check_processed = 0;

cases_list_ai_used = cases_list_ai;
for k=1:size(cases_list_ai, 1)
    a_case = cases_list_ai{k};
    a_case_new = strrep(a_case, '/data/raw_data/retro_cine_3T/raw_selected/', '');
    cases_list_ai_used{k} = a_case_new;
end

num_cases = size(cases_list_ai_used, 1)

script_dir = '/data/raw_data/retro_cine_3T/local_scripts_6/'
mkdir(script_dir)
delete(fullfile(script_dir, '*.sh'))

setenv('model_dir', '/data/raw_data/retro_cine_3T/models/')
setenv('data_dir', '/data/raw_data/retro_cine_3T/raw_selected/')
setenv('output_dir', '/data/raw_data/retro_cine_3T/raw_selected/')
setenv('code_dir', '/home/xueh/mrprogs/imagingfm_BTCHW/')
setenv('batch_size', '1')

check_processed = 1;
added_noise_sd = 0.1
rep = 1
input_model_list = [];
input_model_list = {
    "hrnet_200m", "hrnet_200m.pth";   
    };

run_file_name = fullfile(script_dir, ['run_snr_data_R_test.sh'])
gpu_lists={'0', '1', '2', '3'}
model_list = run_model_inference_v2(run_file_name, 1, cases_list_ai_used, gpu_lists, check_processed, added_noise_sd, rep, input_model_list);

% run to get clean im
script_dir = '/data/raw_data/retro_cine_3T/local_scripts_clean_im2/'
mkdir(script_dir)
delete(fullfile(script_dir, '*.sh'))
model_file = '/data/models/paper_v3/model_search2_hrnet__whitney15_4xG8_MI300X_epoch_162_98p/whitney15_hrnet__deployed_200m_2dt_4xG8/checkpoints/paper-STCNNT_HRNET_T1L1G1T1L1G1T1L1G1T1L1G1T1L1G1T1L1G1T1L1G1T1L1G1_T1L1G1T1L1G1T1L1G1T1L1G1T1L1G1T1L1G1T1L1G1T1L1G1_whitney15_hrnet__deployed_200m_2dt_4xG8_4x192G8-MI300X_NN_32.0_32.0_C-64_complex_residual/checkpoint_epoch_161.pth'
res_dir = 'res_R2'
setenv('data_dir', '/data/raw_data/retro_cine_3T/raw_selected/')
setenv('output_dir', '/data/raw_data/retro_cine_3T/raw_selected/')
setenv('code_dir', '/home/xueh/mrprogs/imagingfm_BTCHW/')
setenv('batch_size', '1')

check_processed = 1;
run_file_name = fullfile(script_dir, ['run_snr_data_clean_im.sh'])
gpu_lists={'0', '1', '2', '3'}

run_a_model_inference(run_file_name, cases_list, gpu_lists, model_file, res_dir, check_processed)

% for tt=1:numel(cases_list)
%     a_case = cases_list{tt}
%     rmdir(fullfile(a_case, '../../res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE'), 's');
% end

% copy for cloud computing

cloud_dir = '/data/raw_data/retro_cine_3T/cloud'
mkdir(cloud_dir)

for tt=1:numel(cases_list_ai)
    a_case = cases_list_ai{tt}
    a_case_dir = a_case(43:end)
    a_case_dst_dir = fullfile(cloud_dir, a_case_dir)

    if ~exist((fullfile(a_case_dst_dir, 'gmap.npy')))
        mkdir(a_case_dst_dir)
    
        disp([num2str(tt) ' - copying to ' a_case_dst_dir])
    
        copyfile(fullfile(a_case, 'input_real.npy'), a_case_dst_dir);
        copyfile(fullfile(a_case, 'input_imag.npy'), a_case_dst_dir);
        copyfile(fullfile(a_case, 'gmap.npy'), a_case_dst_dir);
    end
end

FileList = dir(fullfile(cloud_dir, '**', 'gmap.npy'));

cases_list_ai_cloud = [];
for k=1:numel(FileList)
    cases_list_ai_cloud = [cases_list_ai_cloud; {FileList(k).folder}];
end

cases_list_ai_cloud_used = cases_list_ai_cloud;
for k=1:size(cases_list_ai_cloud, 1)
    a_case = cases_list_ai_cloud{k};
    a_case_new = strrep(a_case, '/data/raw_data/retro_cine_3T/cloud/', '');
    cases_list_ai_cloud_used{k} = a_case_new;
end


project = 'v6'

script_dir = ['/data/raw_data/retro_cine_3T/scripts_' project]
delete(fullfile(script_dir, '*.sh'))
script_name = 'run_snr_data'

num_cases = size(cases_list_ai_cloud_used, 1)
nnodes = 27

setenv('model_dir', '/scratch/amlt_code/models/')
setenv('data_dir', '/hot/projects/ifm/paper2/cloud/')
setenv('output_dir', '/hot/projects/ifm/paper2/cloud/')
setenv('local_output_dir', '/data/raw_data/retro_cine_3T/raw_selected/')
setenv('code_dir', '/scratch/amlt_code')
setenv('batch_size', '1')

mkdir(script_dir)

check_processed = 1;
added_noise_sd = 0.1
rep = 1

input_model_list = {
    "hrnet_28m", "hrnet_28m.pth";
    "hrnet_50m",  "hrnet_50m.pth";     
    "hrnet_100m", "hrnet_100m.pth"; 
    "hrnet_200m", "hrnet_200m.pth";   

    "hrnet_TTT", "hrnet_TTT.pth";
    "hrnet_TTTTTT", "hrnet_TTTTTT.pth";
    "hrnet_Swin", "hrnet_Swin.pth";
    "hrnet_ViT3D", "hrnet_ViT3D.pth";
    "hrnet_ViT3D_large", "hrnet_ViT3D_large.pth";
    "hrnet_C3", "hrnet_C3.pth";
    "hrnet_C3_large", "hrnet_C3_large.pth"
    };

B = linspace(1, num_cases, nnodes+1)
B = round(B)

run_file_name = fullfile(script_dir, [script_name '.sh'])
gpu_lists={'0,1', '2,3', '4,5', '6,7', '8,9', '10,11', '12,13', '14,15'}
model_list = run_model_inference_v2(run_file_name, nnodes, cases_list_ai_cloud_used, gpu_lists, check_processed, added_noise_sd, rep, input_model_list)
    
sub_dirs = {'hrnet_28m', 'hrnet_50m', 'hrnet_100m', 'hrnet_200m'}
sub_dirs = {'hrnet_28m', 'hrnet_50m', 'hrnet_100m', 'hrnet_200m', 'hrnet_TTT', 'hrnet_TTTTTT', 'hrnet_Swin', 'hrnet_ViT3D', 'hrnet_ViT3D_large', 'hrnet_C3', 'hrnet_C3_large'}

% copy cloud res to local

for k=1:numel(cases_list_ai_used)
    case_dir = fullfile('/data/raw_data/retro_cine_3T/cloud_res/cloud/', cases_list_ai_used{k});
    if exist(case_dir)
        try
            if exist(fullfile(case_dir, 'gmap.npy'))
                disp(['delete - ' case_dir])
                delete(fullfile(case_dir, 'gmap.npy'));
                delete(fullfile(case_dir, 'input_real.npy'));
                delete(fullfile(case_dir, 'input_imag.npy'));
            end
        catch
        end
    end
end

for i=1:size(cases_list_ai_used, 1)

    src_dir = fullfile('/data/raw_data/retro_cine_3T/cloud_res/cloud', cases_list_ai_used{i});
    [sub_dirs, num_sub_dirs] = FindSubDirs(fullfile('/data/raw_data/retro_cine_3T/cloud_res/cloud', cases_list_ai_used{i}));

    for s=1:num_sub_dirs

        src_run_dir = fullfile(src_dir, sub_dirs{s});

        if exist(fullfile(src_run_dir, 'output_real.npy')) & exist(fullfile(src_run_dir, 'output_imag.npy'))
            dst_dir = fullfile('/data/raw_data/retro_cine_3T/raw_selected', cases_list_ai_used{i}, sub_dirs{s});
    
            if ~exist(fullfile(dst_dir, 'output_real.npy')) | ~exist(fullfile(dst_dir, 'output_imag.npy'))
                if ~exist(dst_dir)
                    mkdir(dst_dir)
                end
    
                disp(['copy ' cases_list_ai_used{i} ' for ' sub_dirs{s}]);
                movefile(fullfile(src_run_dir, 'output_real.npy'), dst_dir);
                movefile(fullfile(src_run_dir, 'output_imag.npy'), dst_dir);
            else
                delete(fullfile(src_run_dir, 'output_real.npy'));
                delete(fullfile(src_run_dir, 'output_imag.npy'));
            end
        end
    end
end

% remove results
sub_dirs = {'hrnet_Swin', 'hrnet_ViT3D', 'hrnet_ViT3D_large', 'hrnet_C3', 'hrnet_C3_large'}

for k=1:numel(cases_list_ai_used)
    case_dir = fullfile('/data/raw_data/retro_cine_3T/raw_selected/', cases_list_ai_used{k});
    if exist(case_dir)
        for t=1:numel(sub_dirs)
            slc_dir = fullfile(case_dir, sub_dirs{t});
            if exist(fullfile(slc_dir, 'output_real.npy'))
                slc_dir
                a = complex(readNPY(fullfile(slc_dir, 'output_real')), readNPY(fullfile(slc_dir, 'output_imag')));
                writeNPY(abs(a), fullfile(slc_dir, 'output.npy'));
                delete(fullfile(slc_dir, 'output_real.npy'));
                delete(fullfile(slc_dir, 'output_imag.npy'));
            end
        end
    end
end


% compute ssim, psnr etc.

N = size(cases_list_ai, 1)

c_list = []

ssim_28m = []
ssim_50m = []
ssim_100m = []
ssim_200m = []
ssim_TTT = []
ssim_TTTTTT = []
ssim_Swin = []
ssim_ViT3D = []
ssim_ViT3D_large = []
ssim_C3 = []
ssim_C3_large = []

psnr_28m = []
psnr_50m = []
psnr_100m = []
psnr_200m = []
psnr_TTT = []
psnr_TTTTTT = []
psnr_Swin = []
psnr_ViT3D = []
psnr_ViT3D_large = []
psnr_C3 = []
psnr_C3_large = []

sub_dirs = {'hrnet_28m', 'hrnet_50m', 'hrnet_100m', 'hrnet_200m'}
sub_dirs = {'hrnet_28m', 'hrnet_50m', 'hrnet_100m', 'hrnet_200m', 'hrnet_TTT', 'hrnet_TTTTTT', 'hrnet_Swin', 'hrnet_ViT3D', 'hrnet_ViT3D_large', 'hrnet_C3', 'hrnet_C3_large'}

only_metrics = 1;

for tt=1:N
    disp([' ================ ' num2str(tt) ' out of ' num2str(N) ' ================ '])
    a_case = cases_list_ai{tt};
    try
        [ssims, psnrs, noise_sigmas, gmap, clean_im, input, outputs] = compare_models_metric(a_case, sub_dirs, 0, only_metrics);
        
        ssim_28m = [ssim_28m; ssims(1, :)];
        ssim_50m = [ssim_50m; ssims(2, :)];
        ssim_100m = [ssim_100m; ssims(3, :)];

        psnr_28m = [psnr_28m; psnrs(1, :)];
        psnr_50m = [psnr_50m; psnrs(2, :)];
        psnr_100m = [psnr_100m; psnrs(3, :)];

        ssim_200m = [ssim_200m; ssims(4, :)];
        psnr_200m = [psnr_200m; psnrs(4, :)];

        ssim_TTT = [ssim_TTT; ssims(5, :)];
        psnr_TTT = [psnr_TTT; psnrs(5, :)];

        ssim_TTTTTT = [ssim_TTTTTT; ssims(6, :)];
        psnr_TTTTTT = [psnr_TTTTTT; psnrs(6, :)];

        ssim_Swin = [ssim_Swin; ssims(7, :)];
        psnr_Swin = [psnr_Swin; psnrs(7, :)];

        ssim_ViT3D = [ssim_ViT3D; ssims(8, :)];
        psnr_ViT3D = [psnr_ViT3D; psnrs(8, :)];

        ssim_ViT3D_large = [ssim_ViT3D_large; ssims(9, :)];
        psnr_ViT3D_large = [psnr_ViT3D_large; psnrs(9, :)];

        ssim_C3 = [ssim_C3; ssims(10, :)];
        psnr_C3 = [psnr_C3; psnrs(10, :)];

        ssim_C3_large = [ssim_C3_large; ssims(11, :)];
        psnr_C3_large = [psnr_C3_large; psnrs(11, :)];

        c_list = [c_list; {a_case}];
    catch
        disp(['error - ' a_case]);
    end
end


case_res = table(c_list, ssim_28m, ssim_50m, ssim_100m, ssim_200m, ...
    ssim_TTT, ssim_TTTTTT, ssim_Swin, ssim_ViT3D, ssim_ViT3D_large, ssim_C3, ssim_C3_large, ... 
    psnr_28m, psnr_50m, psnr_100m, psnr_200m, ...
    psnr_TTT, psnr_TTTTTT, psnr_Swin, psnr_ViT3D, psnr_ViT3D_large, psnr_C3, psnr_C3_large);

ind = find(snr_levels>0)
x = snr_levels(ind)
size(x)

% IT vs. others
h_fig = figure('WindowState', 'maximized');
y = cat(3, ssim_28m, ssim_50m, ssim_TTT, ssim_TTTTTT, ssim_Swin, ssim_ViT3D, ssim_ViT3D_large, ssim_C3, ssim_C3_large);
y_m = mean(y, 1);
y_m = squeeze(y_m);
y_m = y_m(ind, :);
size(y_m)
subplot(1, 2, 1)
hold on
h = plot(x, y_m, 'LineWidth', 2.0);
h(1).LineWidth = 4; h(1).Color = [0.5 0.5 0];
h(2).LineWidth = 4; h(2).Color = [0.615 0 0];

h(3).Color = [0.0 0.65 0.0];
h(3).LineStyle = '--';
h(4).Color = [0.0 1.0 0.0];
h(4).LineStyle = '--';

h(5).Color = [0.0 0.0 0.0]; h(7).LineStyle = '-';

h(6).Color = [0.0 0.0 0.65]; h(8).LineStyle = ':';
h(7).Color = [0.0 0.0 1.0]; h(9).LineStyle = ':';

h(8).Color = 0.7*[0.929000000000000 0.694000000000000 0.125000000000000]; h(8).LineStyle = '-.';
h(9).Color = [0.929000000000000 0.694000000000000 0.125000000000000]; h(9).LineStyle = '-.';

grid on
box on
hold off
legend('IT-27m', 'IT-55m', 'CNNT-27m', 'CNNT-55m', 'Swin3D-55m', 'ViT3D-27m', 'ViT3D-55m','Conv3D-23m', 'Conv3D-45m', 'Location','southeast')
xscale log
yscale log
xticks(snr_levels)
xlabel('SNR, input')
ylabel('SSIM after model')
fontsize(26,"points")
title('Compare IT models with others, SSIM')

y = cat(3, psnr_28m, psnr_50m, psnr_TTT, psnr_TTTTTT, psnr_Swin, psnr_ViT3D, psnr_ViT3D_large, psnr_C3, psnr_C3_large);
y_m = mean(y, 1);
y_m = squeeze(y_m);
y_m = y_m(ind, :);
size(y_m)

subplot(1, 2, 2)
hold on
h = plot(x, y_m, 'LineWidth', 2.0);
h(1).LineWidth = 4; h(1).Color = [0.5 0.5 0];
h(2).LineWidth = 4; h(2).Color = [0.615 0 0];

h(3).Color = [0.0 0.65 0.0];
h(3).LineStyle = '--';
h(4).Color = [0.0 1.0 0.0];
h(4).LineStyle = '--';

h(5).Color = [0.0 0.0 0.0]; h(7).LineStyle = '-';

h(6).Color = [0.0 0.0 0.65]; h(8).LineStyle = ':';
h(7).Color = [0.0 0.0 1.0]; h(9).LineStyle = ':';

h(8).Color = 0.7*[0.929000000000000 0.694000000000000 0.125000000000000]; h(8).LineStyle = '-.';
h(9).Color = [0.929000000000000 0.694000000000000 0.125000000000000]; h(9).LineStyle = '-.';
grid on
box on
hold off
legend('IT-27m', 'IT-55m', 'CNNT-27m', 'CNNT-55m', 'Swin3D-55m', 'ViT3D-27m', 'ViT3D-55m','Conv3D-23m', 'Conv3D-45m', 'Location','southeast')
xscale log
xticks(snr_levels)
xlabel('SNR, input')
ylabel('PSNR after model')
fontsize(26,"points")
title('Compare IT models with others, PSNR')

saveas(h_fig, fullfile('/data/raw_data/retro_cine_3T/', 'IT_vs_others.fig'))

% ITs
h_fig = figure('WindowState', 'maximized');
y = cat(3, ssim_28m, ssim_50m, ssim_100m, ssim_200m);
y_m = mean(y, 1);
y_m = squeeze(y_m);
y_m = y_m(ind, :);
size(y_m)
subplot(1, 2, 1)
hold on
h = plot(x, y_m, 'LineWidth', 4.0);
h(1).LineWidth = 2; h(1).Color = [0.5 0.5 0];
h(2).LineWidth = 2; h(2).Color = [0.615 0 0];
h(3).LineWidth = 4; h(3).Color = [0.75 0 0];
h(4).LineWidth = 4; h(4).Color = [1.0 0 0];

grid on
box on
hold off
legend('IT-27m', 'IT-55m', 'IT-109m', 'IT-218m', 'Location','southeast')
xscale log
yscale log
xticks(snr_levels)
xlabel('SNR, input')
ylabel('SSIM after model')
fontsize(26,"points")
title('Compare models with different sizes, SSIM')

y = cat(3, psnr_28m, psnr_50m, psnr_100m, psnr_200m);
y_m = mean(y, 1);
y_m = squeeze(y_m);
y_m = y_m(ind, :);
size(y_m)

subplot(1, 2, 2)
hold on
h = plot(x, y_m, 'LineWidth', 4.0);
h(1).LineWidth = 2; h(1).Color = [0.5 0.5 0];
h(2).LineWidth = 2; h(2).Color = [0.615 0 0];
h(3).LineWidth = 4; h(3).Color = [0.75 0 0];
h(4).LineWidth = 4; h(4).Color = [1.0 0 0];
grid on
box on
hold off
legend('IT-27m', 'IT-55m', 'IT-109m', 'IT-218m', 'Location','southeast')
xscale log
yscale log
xticks(snr_levels)
xlabel('SNR, input')
ylabel('PSNR after model')
fontsize(26,"points")
title('Compare models with different sizes, PSNR')

saveas(h_fig, fullfile('/data/raw_data/retro_cine_3T/', 'IT_models.fig'))

save /data/raw_data/retro_cine_3T/case_res.mat

% create the dicom style data
for n=1:fig_num
    [p, fname, ext] = fileparts(fignames{n});
    disp([num2str(n) ' out of ' num2str(fig_num) ' - ' fname]);

    for kk=4:4
        if kk==1
            final_dir = '/data/raw_data/retro_cine_3T/dicom_final/final_28m/';
            %mkdir(final_dir)
            model_str = 'hrnet_28m'
        end
        if kk==2
            final_dir = '/data/raw_data/retro_cine_3T/dicom_final/final_50m/';
            %mkdir(final_dir)
            model_str = 'hrnet_50m'
        end
        if kk==3
            final_dir = '/data/raw_data/retro_cine_3T/dicom_final/final_100m/';
            %mkdir(final_dir)
            model_str = 'hrnet_100m'
        end
        if kk==4
            final_dir = '/data/raw_data/retro_cine_3T/dicom_final/final_200m/';
            %mkdir(final_dir)
            model_str = 'hrnet_200m'
        end    
    
        ind = find(fname=='_');
    
        case_name = fname(1:ind(1)-1)
    
        final_case_dir = fullfile(final_dir, case_name)
        mkdir(final_case_dir)
    
        % if exist(fullfile(final_case_dir, 'gt_sax_imag.npy'))
        %     continue
        % end
    
        d_dir = fullfile(dst_dir, case_name, 'ch2');
        [h5_dirs, num_dirs] = FindSubDirs(d_dir);
        ch2_dir = fullfile(d_dir, h5_dirs{1})
    
        d_dir = fullfile(dst_dir, case_name, 'ch4');
        [h5_dirs, num_dirs] = FindSubDirs(d_dir);
        ch4_dir = fullfile(d_dir, h5_dirs{1})
    
        d_dir = fullfile(dst_dir, case_name, 'sax');
        [h5_dirs, num_dirs] = FindSubDirs(d_dir);
        sax_dir = fullfile(d_dir, h5_dirs{1})
    
        ch2_rec = load(fullfile(ch2_dir, res_name, 'model', 'for_model.mat'));
        ch4_rec = load(fullfile(ch4_dir, res_name, 'model', 'for_model.mat'));
        sax_rec = load(fullfile(sax_dir, res_name, 'model', 'for_model.mat'));
    
        ch2_model_dir = fullfile(ch2_dir, res_name, 'model');
        ch4_model_dir = fullfile(ch4_dir, res_name, 'model');
        sax_model_dir = fullfile(sax_dir, res_name, 'model');
    
        % data_ch2 = complex(readNPY(fullfile(ch2_model_dir, 'input_real.npy')), readNPY(fullfile(ch2_model_dir, 'input_imag.npy')));
        % data_ch4 = complex(readNPY(fullfile(ch4_model_dir, 'input_real.npy')), readNPY(fullfile(ch4_model_dir, 'input_imag.npy')));
        data_sax_raw = complex(readNPY(fullfile(sax_model_dir, 'input_real.npy')), readNPY(fullfile(sax_model_dir, 'input_imag.npy')));

        data_ch2 = complex(readNPY(fullfile(ch2_model_dir, 'res_R2', 'output_real.npy')), readNPY(fullfile(ch2_model_dir, 'res_R2', 'output_imag.npy')));
        data_ch4 = complex(readNPY(fullfile(ch4_model_dir, 'res_R2', 'output_real.npy')), readNPY(fullfile(ch4_model_dir, 'res_R2', 'output_imag.npy')));
        data_sax = complex(readNPY(fullfile(sax_model_dir, 'res_R2', 'output_real.npy')), readNPY(fullfile(sax_model_dir, 'res_R2', 'output_imag.npy')));
    
        gmap_sax = readNPY(fullfile(sax_model_dir, 'gmap'));

        size(data_ch2)
        size(data_ch4)
        size(data_sax)
    
        SLC = size(data_sax, 4);
    
        gt_h_4ch = ch4_rec.header{1};
        gt_h_2ch = ch2_rec.header{1};
        
        gt_h_2ch.matrix_size = size(data_ch2, [2 1 3]);
        gt_h_4ch.matrix_size = size(data_ch4, [2 1 3]);
        
        gt_2ch = squeeze(abs(data_ch2)); size(gt_2ch)
        gt_4ch = squeeze(abs(data_ch4)); size(gt_4ch)
        gt_sax = squeeze(abs(data_sax)); size(gt_sax)
        gt_sax_raw = squeeze(abs(data_sax_raw)); size(gt_sax_raw)
    
        writeNPY(single(gt_2ch), fullfile(final_case_dir, 'gt_2ch.npy'));
        writeNPY(single(gt_4ch), fullfile(final_case_dir, 'gt_4ch.npy'));
        % writeNPY(single(gt_sax), fullfile(final_case_dir, 'gt_sax.npy'));
        % writeNPY(single(gt_sax_raw), fullfile(final_case_dir, 'gt_sax_raw.npy'));
        writeNPY(gmap_sax, fullfile(final_case_dir, 'gt_gmap_sax.npy'));

        h = figure; imagescn(gt_2ch, [0 4*mean(gt_2ch(:))], [], [12], 3);
        saveas(h, fullfile(final_case_dir, ['gt_2ch.jpg']), 'jpg');

        h = figure; imagescn(gt_4ch, [0 4*mean(gt_4ch(:))], [], [12], 3);
        saveas(h, fullfile(final_case_dir, ['gt_4ch.jpg']), 'jpg');

        h = figure; imagescn(gt_sax, [0 4*mean(gt_sax(:))], [], [12], 3);
        saveas(h, fullfile(final_case_dir, ['gt_sax.jpg']), 'jpg');

        c_sax = single(squeeze(data_sax));
        writeNPY(single(real(c_sax)), fullfile(final_case_dir, 'gt_sax_real.npy'));
        writeNPY(single(imag(c_sax)), fullfile(final_case_dir, 'gt_sax_imag.npy'));

        c_sax_raw = single(squeeze(data_sax_raw));
        writeNPY(single(real(c_sax_raw)), fullfile(final_case_dir, 'gt_sax_raw_real.npy'));
        writeNPY(single(imag(c_sax_raw)), fullfile(final_case_dir, 'gt_sax_raw_imag.npy'));
    
        h_gt_4ch = ComputeDicomCoordFromGtOffline(gt_h_4ch.PatientPosition, gt_h_4ch.read_dir, gt_h_4ch.phase_dir, double(gt_h_4ch.FOV), double(gt_h_4ch.matrix_size), size(gt_4ch, 2), size(gt_4ch, 1));
        h_gt_2ch = ComputeDicomCoordFromGtOffline(gt_h_2ch.PatientPosition, gt_h_2ch.read_dir, gt_h_2ch.phase_dir, double(gt_h_2ch.FOV), double(gt_h_2ch.matrix_size), size(gt_2ch, 2), size(gt_2ch, 1));
    
        h_gt_sax = [];
        for slc_sax=1:SLC
            gt_h_sax = sax_rec.header{slc_sax, 1, 1};
            gt_h_sax.matrix_size = size(data_sax, [2 1 3]); gt_h_sax.matrix_size(3) = 1;
        
            
            h_gt_sax{slc_sax} = ComputeDicomCoordFromGtOffline(gt_h_sax.PatientPosition, gt_h_sax.read_dir, gt_h_sax.phase_dir, double(gt_h_sax.FOV), double(gt_h_sax.matrix_size), size(gt_sax, 2), size(gt_sax, 1));
        end

        h = figure('Name', case_name);
        ah = axes;
        plot2DDicomSlice(gt_4ch(:,:,1), h_gt_4ch, ah, -1);
        plot2DDicomSlice(gt_2ch(:,:,1), h_gt_2ch, ah, -1);
    
        gt_sax_picked = squeeze(abs(data_sax(:,:,:,round(SLC/2))));
        plot2DDicomSlice(gt_sax_picked(:,:,1), h_gt_sax{round(SLC/2)}, ah, -1);
    
        saveas(h, fullfile(final_case_dir, [case_name '_dicom_plot.fig']), 'fig');
        saveas(h, fullfile(final_case_dir, [case_name '_dicom_plot.jpg']), 'jpg');
    
        close(h)
    
        % try
            ch2_res_dir = fullfile(ch2_dir, res_name, 'model', 'res', 'slc_0');
            ch4_res_dir = fullfile(ch4_dir, res_name, 'model', 'res', 'slc_0');
            sax_res_dir = fullfile(sax_dir, res_name, 'model', 'res');
        
            ch2_noise_sigmas = readNPY(fullfile(ch2_res_dir, 'noise_sigmas.npy'));
            ch2_input = complex(readNPY(fullfile(ch2_res_dir, 'input_real.npy')), readNPY(fullfile(ch2_res_dir, 'input_imag.npy')));
            ch2_output = complex(readNPY(fullfile(ch2_res_dir, model_str, 'output_real.npy')), readNPY(fullfile(ch2_res_dir, model_str, 'output_imag.npy')));
            size(ch2_input)
            size(ch2_output)
        
            for k=1:numel(ch2_noise_sigmas)-1
                ch2_input(:,:,:,k) = ch2_input(:,:,:,k) * ch2_noise_sigmas(k);
                ch2_output(:,:,:,k) = ch2_output(:,:,:,k) * ch2_noise_sigmas(k);
            end
        
            h = figure;
            imagescn(cat(4, abs(ch2_input), abs(ch2_output)), [0 5*mean(gt_2ch(:))], [2 numel(ch2_noise_sigmas)-1], [32], 3);
            saveas(h, fullfile(final_case_dir, ['ch2_input_output.jpg']), 'jpg');
        
            ch4_noise_sigmas = readNPY(fullfile(ch4_res_dir, 'noise_sigmas.npy'));
            ch4_input = complex(readNPY(fullfile(ch4_res_dir, 'input_real.npy')), readNPY(fullfile(ch4_res_dir, 'input_imag.npy')));
            ch4_output = complex(readNPY(fullfile(ch4_res_dir, model_str, 'output_real.npy')), readNPY(fullfile(ch4_res_dir, model_str, 'output_imag.npy')));
            size(ch4_input)
            size(ch4_output)
        
            for k=1:numel(ch4_noise_sigmas)-1
                ch4_input(:,:,:,k) = ch4_input(:,:,:,k) * ch4_noise_sigmas(k);
                ch4_output(:,:,:,k) = ch4_output(:,:,:,k) * ch4_noise_sigmas(k);
            end
        
            h = figure;
            imagescn(cat(4, abs(ch4_input), abs(ch4_output)), [0 5*mean(gt_4ch(:))], [2 numel(ch4_noise_sigmas)-1], [32], 3);
            saveas(h, fullfile(final_case_dir, ['ch4_input_output.jpg']), 'jpg');
        
            sax_noise_sigmas = [];
            sax_input = [];
            sax_output = [];
            sax_gmaps_input = [];
            for slc_sax=1:SLC
                disp(['load sax slice ' num2str(slc_sax) ' ... ']);
                slc_dir = fullfile(sax_res_dir, ['slc_' num2str(slc_sax-1)]);
                sax_noise_sigmas(:, slc_sax) = readNPY(fullfile(slc_dir, 'noise_sigmas.npy'));
        
                sax_input(:,:,:,:,slc_sax) = complex(readNPY(fullfile(slc_dir, 'input_real.npy')), readNPY(fullfile(slc_dir, 'input_imag.npy')));
                sax_output(:,:,:,:,slc_sax) = complex(readNPY(fullfile(slc_dir, model_str, 'output_real.npy')), readNPY(fullfile(slc_dir, model_str, 'output_imag.npy')));

                a_gmap = readNPY(fullfile(slc_dir, 'gmap'));
                sax_gmaps_input(:,:,slc_sax) = a_gmap(:,:,1);
            end
        
            sax_input = permute(sax_input, [1 2 3 5 4]);
            sax_output = permute(sax_output, [1 2 3 5 4]);
        
            size(sax_input)
            size(sax_output)
        
            for slc_sax=1:SLC
                for k=1:size(sax_noise_sigmas,1)-1
                    r = sax_noise_sigmas(k, slc_sax);
                    sax_input(:,:,:,slc_sax, k) = sax_input(:,:,:,slc_sax, k) * r;
                    sax_output(:,:,:,slc_sax, k) = sax_output(:,:,:,slc_sax, k) * r;
                end
            end
        
            % sax_input = readNPY('sax_input.npy'); sax_output = readNPY('sax_output.npy'); SLC = size(sax_input, 4); 
    
            slc_picked = round(SLC/2)
            h = figure;
            sax_im = cat(4, abs(sax_input(:,:,:,slc_picked,:)), abs(sax_output(:,:,:,slc_picked,:)));
            sax_im2 = reshape(sax_im, [size(sax_im,1:3) size(sax_im, 4)*size(sax_im,5)]);
            size(sax_im)
            imagescn(squeeze(sax_im), [0 5*mean(gt_sax(:))], [2 size(sax_noise_sigmas,1)-1], [32], 3);
            saveas(h, fullfile(final_case_dir, ['sax_input_output.jpg']), 'jpg');
        
            N_fig = size(sax_noise_sigmas,1)-1;
        
            h = figure;    
            imagescn(sax_im2(:,:,:, end:-1:1), [0 5*mean(gt_sax(:))], [5 8], [16], 3);
            saveas(h, fullfile(final_case_dir, ['sax_input_output_rows.jpg']), 'jpg');
        
            save(fullfile(final_case_dir, 'dicom_headers.mat'), 'h_gt_4ch', 'h_gt_2ch', 'h_gt_sax');
            % writeNPY(single(gt_2ch), fullfile(final_case_dir, 'gt_2ch.npy'));
            % writeNPY(single(gt_4ch), fullfile(final_case_dir, 'gt_4ch.npy'));
            % writeNPY(single(gt_sax), fullfile(final_case_dir, 'gt_sax.npy'));
        
            % writeNPY(single(real(ch2_input)), fullfile(final_case_dir, 'ch2_input_real.npy'));
            % writeNPY(single(imag(ch2_input)), fullfile(final_case_dir, 'ch2_input_imag.npy'));
            % 
            % writeNPY(single(real(ch4_input)), fullfile(final_case_dir, 'ch4_input_real.npy'));
            % writeNPY(single(imag(ch4_input)), fullfile(final_case_dir, 'ch4_input_imag.npy'));
            % 
            % writeNPY(single(real(sax_input)), fullfile(final_case_dir, 'sax_input_real.npy'));
            % writeNPY(single(imag(sax_input)), fullfile(final_case_dir, 'sax_input_imag.npy'));
            % 
            % writeNPY(single(real(ch2_output)), fullfile(final_case_dir, 'ch2_output_real.npy'));
            % writeNPY(single(imag(ch2_output)), fullfile(final_case_dir, 'ch2_output_imag.npy'));
            % 
            % writeNPY(single(real(ch4_output)), fullfile(final_case_dir, 'ch4_output_real.npy'));
            % writeNPY(single(imag(ch4_output)), fullfile(final_case_dir, 'ch4_output_imag.npy'));

            writeNPY(single(real(sax_output)), fullfile(final_case_dir, 'sax_output_real.npy'));
            writeNPY(single(imag(sax_output)), fullfile(final_case_dir, 'sax_output_imag.npy'));
        
            writeNPY(single(abs(ch2_input)), fullfile(final_case_dir, 'ch2_input.npy'));   
            writeNPY(single(abs(ch4_input)), fullfile(final_case_dir, 'ch4_input.npy'));
            writeNPY(single(abs(sax_input)), fullfile(final_case_dir, 'sax_input.npy'));
        
            writeNPY(single(abs(ch2_output)), fullfile(final_case_dir, 'ch2_output.npy'));  
            writeNPY(single(abs(ch4_output)), fullfile(final_case_dir, 'ch4_output.npy'));
            %writeNPY(single(abs(sax_output)), fullfile(final_case_dir, 'sax_output.npy'));

            writeNPY(sax_gmaps_input, fullfile(final_case_dir, 'sax_gmaps_input.npy'));
    
            writeNPY(single(sax_noise_sigmas), fullfile(final_case_dir, 'sax_noise_sigmas.npy'));
            writeNPY(single(ch4_noise_sigmas), fullfile(final_case_dir, 'ch4_noise_sigmas.npy'));
            writeNPY(single(ch2_noise_sigmas), fullfile(final_case_dir, 'ch2_noise_sigmas.npy'));
        
            closeall
        % catch
        %     disp(['failed for ' final_case_dir])
        %     continue
        % end
    end
end

%% cine ai contouring results

% get the masks and computing CNR
raw_28m_dir = '/data/raw_data/retro_cine_3T/dicom_final/final_28m/'
res_28m_dir = '/data/raw_data/retro_cine_3T/dicom_final/Hui_denoise_28M_7mar/'

raw_50m_dir = '/data/raw_data/retro_cine_3T/dicom_final/final_50m/'
res_50m_dir = '/data/raw_data/retro_cine_3T/dicom_final/results_50m/'

raw_100m_dir = '/data/raw_data/retro_cine_3T/dicom_final/final_100m/'
res_100m_dir = '/data/raw_data/retro_cine_3T/dicom_final/Hui_denoise_100M_14mar/'

raw_200m_dir = '/data/raw_data/retro_cine_3T/dicom_final/final_200m/'
res_200m_dir = '/data/raw_data/retro_cine_3T/dicom_final/Hui_denoise_200M_7mar/'

[cases_dirs, num_dirs] = FindSubDirs(raw_200m_dir);
num_dirs

%% for snr, cnr and signal levels
[cases_dirs, num_dirs] = FindSubDirs(res_200m_dir);
[cases_dirs, num_dirs] = FindSubDirs(res_28m_dir);
[cases_dirs, num_dirs] = FindSubDirs(res_100m_dir);


num_dirs

res_dir = res_200m_dir
for i=1:num_dirs
    case_dir = fullfile(res_dir, cases_dirs{i})
    [names, num_gz] = findFILE(case_dir, '*.npy.gz');
    for nn=1:num_gz
        gunzip(names{nn}, case_dir);
    end
end


for i=1:num_dirs
    case_dir = fullfile(res_dir, cases_dirs{i})
    copyfile(fullfile(case_dir, 'SNR_GT__rounded_masks.npy'), fullfile('/data/raw_data/retro_cine_3T/dicom_final/final_200m', cases_dirs{i}))
end


% signal levels
num_levels = numel(snr_levels)
c_list = []

gt_bp_signals = [];
gt_myo_signals = [];

bp_signals_input = [];
myo_signals_input = [];

model_bp_signals = [];
model_myo_signals = [];

sax_noise_sigmas_list = [];

input_cnrs_slices = [];
gt_cnrs_slices = [];
model_cnrs_slices = [];
gt_bp_signals_slices = [];
gt_myo_signals_slices = [];
model_bp_signals_slices = [];
model_myo_signals_slices = [];

cnr_input_all = cell(num_levels, 1);
cnr_gt_all = cell(num_levels, 1);
cnr_model_all = cell(num_levels, 1);

bp_signals_input_all = cell(num_levels, 1);
bp_signals_gt_all = cell(num_levels, 1);
bp_signals_model_all = cell(num_levels, 1);

myo_signals_input_all = cell(num_levels, 1);
myo_signals_gt_all = cell(num_levels, 1);
myo_signals_model_all = cell(num_levels, 1);

for i=1:num_dirs

    disp(['--> ' num2str(i) '  out of ' num2str(num_dirs) ' - ' cases_dirs{i}])

    case_dir = fullfile(raw_200m_dir, cases_dirs{i})
    seg_masks = readNPY(fullfile(case_dir, 'SNR_GT__rounded_masks.npy'));
    seg_masks = permute(seg_masks, [3 4 2 1]);

    %gt_sax = readNPY(fullfile(case_dir, 'gt_sax.npy'));
    sax_input = readNPY(fullfile(case_dir, 'sax_input.npy'));
    %sax_output = readNPY(fullfile(case_dir, 'sax_output.npy'));

    SLC = size(sax_input, 4);

    rec_file = fullfile(case_dir, 'record.mat');
    if exist(rec_file)
        disp(['loading ' rec_file ' ... ']);
        load(rec_file)

        sax_noise_sigmas_list = [sax_noise_sigmas_list; {sax_noise_sigmas}];
        c_list = [c_list; str2num(cases_dirs{i})];

        gt_bp_signals = [gt_bp_signals; bp_signal];
        gt_myo_signals = [gt_myo_signals; myo_signal];

        model_bp_signals = [model_bp_signals; bp_signal_res(:)'];
        model_myo_signals = [model_myo_signals; myo_signal_res(:)'];

        bp_signals_input = [bp_signals_input; bp_signal_input'];
        myo_signals_input = [myo_signals_input; myo_signal_input'];

        input_cnrs_slices = [input_cnrs_slices; {cnr_input_slices}];
        gt_cnrs_slices = [gt_cnrs_slices; {cnr_gt_slices}];
        model_cnrs_slices = [model_cnrs_slices; {cnr_output_slices}];
        gt_bp_signals_slices = [gt_bp_signals_slices; {bp_signals_slices_gt}];
        gt_myo_signals_slices = [gt_myo_signals_slices; {myo_signals_slices_gt}];
        model_bp_signals_slices = [model_bp_signals_slices; {bp_signals_slices_res}];
        model_myo_signals_slices = [model_myo_signals_slices; {myo_signals_slices_res}];

        slc_ind = find(cnr_gt_slices(1,:)>0)
        for n=1:num_levels
            cnr_input_all{n} = [cnr_input_all{n} cnr_input_slices(n, slc_ind)];
            cnr_gt_all{n} = [cnr_gt_all{n} cnr_gt_slices(n, slc_ind)];
            cnr_model_all{n} = [cnr_model_all{n} cnr_output_slices(n, slc_ind)];

            bp_signals_gt_all{n} = [bp_signals_gt_all{n} bp_signals_slices_gt(n, slc_ind)];
            bp_signals_model_all{n} = [bp_signals_model_all{n} bp_signals_slices_res(n, slc_ind)];

            myo_signals_gt_all{n} = [myo_signals_gt_all{n} myo_signals_slices_gt(n, slc_ind)];
            myo_signals_model_all{n} = [myo_signals_model_all{n} myo_signals_slices_res(n, slc_ind)];
        end

        continue
    end

    gt_sax_c = complex(readNPY(fullfile(case_dir, 'gt_sax_real.npy')), readNPY(fullfile(case_dir, 'gt_sax_imag.npy')));
    gt_sax_raw_c = complex(readNPY(fullfile(case_dir, 'gt_sax_raw_real.npy')), readNPY(fullfile(case_dir, 'gt_sax_raw_imag.npy')));

    gt_sax = abs(gt_sax_c);
    gt_sax_raw = abs(gt_sax_raw_c);

    sax_output_c = complex(readNPY(fullfile(case_dir, 'sax_output_real.npy')), readNPY(fullfile(case_dir, 'sax_output_imag.npy')));
    sax_output = abs(sax_output_c);

    gt_gmap_sax = readNPY(fullfile(case_dir, 'gt_gmap_sax.npy'));
    size(gt_gmap_sax)
    %median(gt_gmap_sax(:))

    sax_gmaps_input = readNPY(fullfile(case_dir, 'sax_gmaps_input.npy'));
    size(sax_gmaps_input)
    %median(sax_gmaps_input(:))

    sax_noise_sigmas = readNPY(fullfile(case_dir, 'sax_noise_sigmas.npy'));
    size(sax_noise_sigmas)

    h = load(fullfile(case_dir, 'dicom_headers.mat'));
    n_vec = cross(h.h_gt_sax{1}.ImageOrientationPatient(1:3), h.h_gt_sax{1}.ImageOrientationPatient(4:6));
    slc_locs = zeros(SLC, 1);
    for slc=1:SLC
        slc_locs(slc) = dot(n_vec, h.h_gt_sax{slc}.ImagePositionPatient);
    end
    [slc_locs_sorted, ind_slc] = sort(slc_locs)

    gt_sax = gt_sax(:,:,:,ind_slc);
    gt_sax_c = gt_sax_c(:,:,:,ind_slc);
    gt_sax_raw_c = gt_sax_raw_c(:,:,:,ind_slc);
    sax_input = sax_input(:,:,:,ind_slc,:);
    sax_output = sax_output(:,:,:,ind_slc,:);
    sax_output_c = sax_output_c(:,:,:,ind_slc,:);
    gt_gmap_sax = gt_gmap_sax(:,:,ind_slc);
    sax_gmaps_input = sax_gmaps_input(:,:,ind_slc);
    sax_noise_sigmas = sax_noise_sigmas(:, ind_slc);

    % now the images are from apical to basal, check seg
    seg_endo_area = [];
    for slc=1:SLC
        total_area = -1;
        phs_count = 0;
        for phs=1:size(seg_masks, 3)
            a_seg = seg_masks(:,:,phs, slc);
            endo_phs = find(a_seg==1);
            if ~isempty(endo_phs)
                phs_count = phs_count + 1;
                if (total_area<0)
                    total_area = numel(endo_phs);
                else
                    total_area = total_area + numel(endo_phs);
                end
            end
        end

        if phs_count > 0
            total_area = total_area / phs_count;
            seg_endo_area = [seg_endo_area; total_area];
        end
    end

    if isempty(seg_endo_area)
        continue
    end

    if seg_endo_area(1) > seg_endo_area(end)
        seg_masks = seg_masks(:,:,:,end:-1:1);
    end

    endo_all = find(seg_masks==1);
    epi_all = find(seg_masks==2);

    RO = size(sax_input, 1)
    E1 = size(sax_input, 2)

    Ims = zeros(RO, E1, 3, SLC);
    for slc=1:SLC
        a_seg_masks = seg_masks(:,:,10,slc);
        endo = find(a_seg_masks==1);
        epi = find(a_seg_masks==2);

        im_a = gt_sax(:,:,10,slc);
        im_b = sax_input(:,:,10,slc,end);
        im_c = sax_output(:,:,10,slc,end);        

        if ~isempty(endo)
            im_a(endo) = 32;
            im_b(endo) = 32;
            im_c(endo) = 32;
        end
        if ~isempty(epi)
            im_a(epi) = 128;
            im_b(epi) = 128;
            im_c(epi) = 128;
        end

        Ims(:,:,1,slc) = im_a;
        Ims(:,:,2,slc) = im_b;
        Ims(:,:,3,slc) = im_c;
    end

    h = figure; imagescn(Ims, [0 1.25*mean(sax_input(:))], [ceil(SLC/2) 6], [12]);
    saveas(h, fullfile(case_dir, '/../../seg_mask_plots/', [cases_dirs{i} '_seg_mask_check.jpg']), 'jpg');
    close(h)
    % 
    % continue

    bp_signal_raw = mean(gt_sax_raw(endo_all(:)));
    myo_signal_raw = mean(gt_sax_raw(epi_all(:)));

    bp_signal = mean(gt_sax(endo_all(:)));
    myo_signal = mean(gt_sax(epi_all(:)));    

    bp_signal_res = zeros(num_levels, 1);
    myo_signal_res = zeros(num_levels, 1);

    bp_signal_input = zeros(num_levels, 1);
    myo_signal_input = zeros(num_levels, 1);

    cnr_input_slices = zeros(num_levels, SLC);
    cnr_gt_slices = zeros(num_levels, SLC);
    cnr_output_slices = zeros(num_levels, SLC);
    
    bp_signals_slices_gt = zeros(num_levels, SLC);
    myo_signals_slices_gt = zeros(num_levels, SLC);

    bp_signals_slices_res = zeros(num_levels, SLC);
    myo_signals_slices_res = zeros(num_levels, SLC);

    for n=1:num_levels

        a_sax_input = sax_input(:,:,:,:,n);
        a_sax_output = sax_output(:,:,:,:,n);

        signal_bloods_input = [];
        signal_myos_input = [];

        for slc=1:SLC
            nn_sigma = sax_noise_sigmas(n, slc);
            %a_sax_input = a_sax_input / nn_sigma;
            a_seg_masks = seg_masks(:,:,:,slc);
            endo = find(a_seg_masks==1);
            epi = find(a_seg_masks==2);

            a_sax_input_slc = a_sax_input(:,:,:,slc);
            a_sax_output_slc = a_sax_output(:,:,:,slc);
            gt_sax_slc = gt_sax(:,:,:,slc);

            diff = gt_sax_slc - a_sax_output_slc;

            % figure; imagescn(cat(4, a_sax_input_slc, a_sax_output_slc, a_seg_masks, diff), [], [], [], 3);

            if ~isempty(endo)
                signal_bloods_input = [signal_bloods_input; a_sax_input_slc(endo(:))];
                sd_output_blood = std(diff(endo(:)));
            end
            if ~isempty(epi)
                signal_myos_input = [signal_myos_input; a_sax_input_slc(epi(:))];
                sd_output_myo = std(diff(epi(:)));
            end
            if ~isempty(endo) & ~isempty(epi)
                s_blood = mean(a_sax_input_slc(endo(:)));
                s_myo = mean(a_sax_input_slc(epi(:)));
                cnr_input_slices(n, slc) = (s_blood-s_myo)/nn_sigma;

                s_blood = mean(a_sax_output_slc(endo(:)));
                s_myo = mean(a_sax_output_slc(epi(:)));
                cnr_output_slices(n, slc) = s_blood-s_myo;
                bp_signals_slices_res(n, slc) = s_blood;
                myo_signals_slices_res(n, slc) = s_myo;

                s_blood = mean(gt_sax_slc(endo(:)));
                s_myo = mean(gt_sax_slc(epi(:)));
                bp_signals_slices_gt(n, slc) = s_blood;
                myo_signals_slices_gt(n, slc) = s_myo;
                cnr_gt_slices(n, slc) = s_blood-s_myo;
            end
        end

        bp_signal_input(n) = mean(signal_bloods_input);
        myo_signal_input(n) = mean(signal_myos_input);

        bp_signal_res(n) = mean(abs(a_sax_output(endo_all(:))));
        myo_signal_res(n) = mean(abs(a_sax_output(epi_all(:))));
    end

    sax_noise_sigmas_list = [sax_noise_sigmas_list; {sax_noise_sigmas}];
    c_list = [c_list; str2num(cases_dirs{i})];

    gt_bp_signals = [gt_bp_signals; bp_signal];
    gt_myo_signals = [gt_myo_signals; myo_signal];

    model_bp_signals = [model_bp_signals; bp_signal_res(:)'];
    model_myo_signals = [model_myo_signals; myo_signal_res(:)'];

    bp_signals_input = [bp_signals_input; bp_signal_input'];
    myo_signals_input = [myo_signals_input; myo_signal_input'];

    input_cnrs_slices = [input_cnrs_slices; {cnr_input_slices}];
    gt_cnrs_slices = [gt_cnrs_slices; {cnr_gt_slices}];
    model_cnrs_slices = [model_cnrs_slices; {cnr_output_slices}];
    gt_bp_signals_slices = [gt_bp_signals_slices; {bp_signals_slices_gt}];
    gt_myo_signals_slices = [gt_myo_signals_slices; {myo_signals_slices_gt}];
    model_bp_signals_slices = [model_bp_signals_slices; {bp_signals_slices_res}];
    model_myo_signals_slices = [model_myo_signals_slices; {myo_signals_slices_res}];

    slc_ind = find(cnr_gt_slices(1,:)>0)
    for n=1:num_levels
        cnr_input_all{n} = [cnr_input_all{n} cnr_input_slices(n, slc_ind)];
        cnr_gt_all{n} = [cnr_gt_all{n} cnr_gt_slices(n, slc_ind)];
        cnr_model_all{n} = [cnr_model_all{n} cnr_output_slices(n, slc_ind)];

        bp_signals_gt_all{n} = [bp_signals_gt_all{n} bp_signals_slices_gt(n, slc_ind)];
        bp_signals_model_all{n} = [bp_signals_model_all{n} bp_signals_slices_res(n, slc_ind)];

        myo_signals_gt_all{n} = [myo_signals_gt_all{n} myo_signals_slices_gt(n, slc_ind)];
        myo_signals_model_all{n} = [myo_signals_model_all{n} myo_signals_slices_res(n, slc_ind)];
    end
    save(fullfile(case_dir, 'record.mat'), 'slc_ind', 'sax_noise_sigmas', 'bp_signal', 'myo_signal','bp_signal_res', 'myo_signal_res', 'bp_signal_input', 'myo_signal_input', 'cnr_input_slices', 'cnr_gt_slices', 'cnr_output_slices', 'bp_signals_slices_gt','myo_signals_slices_gt','bp_signals_slices_res','myo_signals_slices_res')
end


for ss=1:3
    %closeall
    
    % ---------------------------------
    if ss ==1
        model_v = cnr_model_all;
        GT_v = cnr_gt_all{1};
        
        y_label_str = "GT-Model, CNR"
        x_label_str = "(GT+Model)/2, CNR"
        h_fig = figure("Name", "CNR, GT vs. model")
        y_lim = [-25 25]
        x_lim = [min(GT_v)*0.95 max(GT_v)*1.1]
        suffix = ''
        
        fig_name = fullfile('/data/raw_data/retro_cine_3T/dicom_final/', get(h_fig, 'Name'))
        h_fig.WindowState = 'maximized';
    end

    % ---------------------------------
    if ss==2
        model_v = bp_signals_model_all;
        GT_v = bp_signals_gt_all{1};
        
        y_label_str = "GT-Model, Signal, blood pool"
        x_label_str = "(GT+Model)/2, Signal, blood pool"
        h_fig = figure("Name", "Signal, blood pool, GT vs. model")
        y_lim = [-25 25]
        x_lim = [min(GT_v)*0.95 max(GT_v)*1.1]
        suffix = ''
        
        fig_name = fullfile('/data/raw_data/retro_cine_3T/dicom_final/', get(h_fig, 'Name'))
        h_fig.WindowState = 'maximized';
    end
    % ---------------------------------
    if ss==3
        model_v = myo_signals_model_all;
        GT_v = myo_signals_gt_all{1};
        
        y_label_str = "GT-Model, Signal, myocardium"
        x_label_str = "(GT+Model)/2, Signal,myocardium"
        h_fig = figure("Name", "Signal, myocardium, GT vs. model")
        y_lim = [-25 25]
        x_lim = [min(GT_v)*0.95 max(GT_v)*1.1]
        suffix = ''
        
        fig_name = fullfile('/data/raw_data/retro_cine_3T/dicom_final/', get(h_fig, 'Name'))
        h_fig.WindowState = 'maximized';
    end
    % ---------------------------------
    
    abs_mean_diff = zeros(num_levels, 1);
    CRs = zeros(num_levels, 2);
    
    symbol = '+'
    markersize = 12
    CR_range = 1.65
    
    for k=num_levels:-1:1
    
        v = model_v{k};
    
        snr_ind = k;
        snr_levels(snr_ind)
    
        a_v = v;
        ind = find(a_v>0.01 & a_v<1e3);
    
        B = a_v(ind);
        A = GT_v(ind);
    
        subplot(2, 5, num_levels-k+1)
        hold on
        [means,diffs,meanDiff,CR,linFit] = BlandAltman2(A', B', 2, symbol, 0, markersize, CR_range)        
        hold off        
        ylabel(y_label_str)
        xlabel(x_label_str)
        box on
        grid on
        ylim(y_lim)
        xlim(x_lim)
        ha = gca;
        ha.YAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
        ha.YAxis.MinorTickValues = y_lim(1):y_lim(2);
        ha.FontSize = 16
    
        [h, p] = ttest(A', B');
        %title(['SNR = ' num2str(snr_levels(17-k+5+1)) ', mean diff = ' num2str(meanDiff, 4) '%, p=' num2str(p, 5)])    
    
        CRs(k, :) = CR;
        abs_mean_diff(snr_ind) = mean(abs(A-B));
    
        title({['SNR = ' num2str(snr_levels(k)) ', mean diff = ' num2str(meanDiff, 4) suffix]; ['90% CR = ' num2str(CR(1), 4) ' to ' num2str(CR(2), 4) suffix]})
    end
    
    fig_name_saved = strrep(strrep(strrep(fig_name, " ", '_'), ',', '_'), '.', '_')
    save(strjoin([fig_name_saved '_CRs.mat'], ''), 'CRs')
    saveas(h_fig, fig_name_saved, 'fig');

    close(h_fig)
end

% plot the cnr example
cd /data/raw_data/retro_cine_3T/dicom_final/final_200m/00048196/
cd /data/raw_data/retro_cine_3T/dicom_final/final_200m/00048233/
cd /data/raw_data/retro_cine_3T/dicom_final/final_200m/00048190/
cd /data/raw_data/retro_cine_3T/dicom_final/final_200m/00048153/

record = load('record')
record.cnr_gt_slices
record.cnr_output_slices

sax_input = readNPY('sax_input.npy');
size(sax_input)

sax_output_c = complex(readNPY('sax_output_real.npy'), readNPY('sax_output_imag.npy'));
sax_output = abs(sax_output_c);

gt_sax_c = complex(readNPY('gt_sax_real.npy'), readNPY('gt_sax_imag.npy'));
gt_sax = abs(gt_sax_c);

seg_masks = readNPY('SNR_GT__rounded_masks.npy');
seg_masks = permute(seg_masks, [3 4 2 1]);
size(seg_masks)



slc = 4

a = record.cnr_gt_slices(:, slc)
b = record.cnr_output_slices(:, slc)
c = record.cnr_input_slices(:, slc)
ns = record.sax_noise_sigmas(:, slc)

c = a./sqrt(1 + ns(1:10).^2)

h = figure('Name', 'input SNR vs. CNR')
ha = gca;
hold on
plot(snr_levels, a, 'k--', 'LineWidth', 4.0)
plot(snr_levels, b, 'r-', 'LineWidth', 2.0)
plot(snr_levels, c, 'b-', 'LineWidth', 2.0)
hold off
grid on
box on
xlim([0 snr_levels(end)])
xlabel('input SNR')
ylabel('CNR')
ha.FontSize = 16
xticks(snr_levels)
xscale log
legend('ground-truch', 'model', 'input')
saveas(h, fullfile('/data/raw_data/retro_cine_3T/dicom_final/CNR_plot.fig'), 'fig')


h = figure('Name', 'CNR, before and after model')
ha = gca;
hold on
plot(c, a, 'k--', 'LineWidth', 4.0)
plot(c, b, 'rx', 'LineWidth', 2.0)
hold off
grid on
box on
xlabel('input CNR')
ylabel('model CNR')
%ylim([80 115])
ha.FontSize = 16
legend('input vs. ground-truch', 'input vs. model')
saveas(h, fullfile('/data/raw_data/retro_cine_3T/dicom_final/CNR_plot_input_vs_model.fig'), 'fig')



im = cat(4,sax_input(:,:,:,slc,:),sax_output(:,:,:,slc,:));
snr_levels([1:4 8])
h=figure; imagescn(im(:,:,:,:,[1:4 8]), [0 8*mean(gt_sax(:))], [2 5], [12], 3)
saveas(h, fullfile('/data/raw_data/retro_cine_3T/dicom_final/CNR_plot_input_model.fig'), 'fig')

h = figure;imagescn(gt_sax(:,:,:,slc), [0 8*mean(gt_sax(:))], [], [12], 3)
saveas(h, fullfile('/data/raw_data/retro_cine_3T/dicom_final/CNR_plot_gt.fig'), 'fig')

a_seg = seg_masks(:,:,:,[end:-1:1]);
a_seg  = a_seg(:,:,:,slc);

endo = find(a_seg==1);
epi = find(a_seg==2);

a_gt_sax = gt_sax(:,:,:,slc);
a_gt_sax(endo) = 32;
a_gt_sax(epi) = 256;
h = figure;imagescn(a_gt_sax, [0 8*mean(gt_sax(:))], [], [12], 3)
saveas(h, fullfile('/data/raw_data/retro_cine_3T/dicom_final/CNR_plot_gt_with_seg.fig'), 'fig')

%% prepare for Expert reading
re_run_flag = 1

ssims = case_res.ssim_200m(:, 4)
ind = find(ssims>0.8);
ind_p = randperm(numel(ind));

ind_picked= ind(ind_p)

cases_experts = case_res.c_list(ind_picked);

expert_dir = '/data/raw_data/retro_cine_3T/dicom_final/expert_cases_200m'
mkdir(expert_dir)

base_name = 'hrnet_200m'
for t=1:size(cases_experts, 1)

    case_dir = cases_experts{t, 1};
    disp([num2str(t) ' - ' case_dir])

    inds = find(case_dir=='/')

    case_name = case_dir(inds(5)+1:inds(6)-1)
    view_str = case_dir(inds(6)+1:inds(7)-1)
    slc_str = case_dir(inds(end)+1:end)

    slc = str2num(slc_str(5:end))

    case_str = [case_name '_' view_str '_' slc_str];

    fig_str = fullfile(expert_dir, case_str)
    fig_name = fullfile(expert_dir, [case_str '.jpg'])
    if exist(fig_name) & ~re_run_flag
        continue
    end

    res_dir = fullfile(case_dir, base_name);

    tic
    noise_sigma = readNPY(fullfile(case_dir, 'noise_sigmas.npy'));

    clean_im = readNPY(fullfile(case_dir, '../../input_real.npy')) + j * readNPY(fullfile(case_dir, '../../input_imag.npy')) ;
    noisy_im = readNPY(fullfile(case_dir, 'input_real.npy')) + j * readNPY(fullfile(case_dir, 'input_imag.npy')) ;

    res = readNPY(fullfile(res_dir, 'output_real.npy')) + j * readNPY(fullfile(res_dir, 'output_imag.npy')) ;
    toc

    tic

    phs = size(res, 3);
    num_snrs = numel(noise_sigma)-1;
    
    r1 = res;
    input1 = noisy_im;
    if ~exist(fig_name) | re_run_flag

        a_clean_im = clean_im(:,:,:,slc+1);

        for k=1:num_snrs
            r1(:,:,:,k) = res(:,:,:,k) * noise_sigma(k);
            input1(:,:,:,k) = noisy_im(:,:,:,k) * noise_sigma(k);
        end

        y = abs(cat(4, input1(:,:,:,end:-1:1), r1(:,:,:,end:-1:1)));
        a = abs(a_clean_im);
        h = figure('name', case_str); imagescn(y, [0 3*mean(abs(a(:)))], [4 size(y, 4)/4], [16], 3);
        %tic; saveas(h, [fig_str '_input_output.fig'], 'fig'); toc
        tic; saveas(h, [fig_str '_input_output.jpg'], 'jpg'); toc
        close(h)       

        y = abs(cat(4, a_clean_im, r1(:,:,:,end:-1:1)));
        a = abs(a_clean_im);
        h = figure('name', case_str); imagescn(y, [0 3*mean(abs(a(:)))], [2 ceil(size(y, 4)/2)], [32], 3);

        axes = findall(h,'type','axes');

        dx = 15
        dy = size(a_clean_im, 1) - 10

        text(axes(11), x, dy, 'GT', 'FontSize', 20, 'Color', [1 1 0], 'FontWeight', 'bold', 'BackgroundColor', [0 0.5 0.5]);
        for k=1:10
            text(axes(k), x, dy, num2str(10-k+1), 'FontSize', 20, 'Color', [1 1 0], 'FontWeight', 'bold', 'BackgroundColor', [0 0.5 0.5]);
        end

        tic; saveas(h, [fig_str '.fig'], 'fig'); toc
        tic; saveas(h, [fig_str '.jpg'], 'jpg'); toc
        close(h)       
    end
end

[case_figs, num_cases] = findFILE(expert_dir, '*.fig');

expert_review_dir = '/data/raw_data/retro_cine_3T/dicom_final/expert_cases_200m_for_review'
mkdir(expert_review_dir)

case_ind = randperm(num_cases)

for t=1:num_cases

    [fpath, case_str, ext] = fileparts(case_figs{case_ind(t), 1});
    disp([num2str(t) ' - ' case_str])

    fig_str = fullfile(expert_dir, case_str)

    fig_name =[fig_str '.fig'];
    h = openfig(fig_name)
    MR_toolbox
    reply = input('Do you want keep this one? y/n [y]:','s');
   if isempty(reply)
      reply = 'y';
   end
   if reply == 'y'
       fig_name
       movefile(fig_name, expert_review_dir);
   end
   closeall
end

[case_figs_selected, num_cases_selected] = findFILE(expert_review_dir, '*.fig');
num_cases_selected

expert_review_gif_dir = '/data/raw_data/retro_cine_3T/dicom_final/expert_cases_200m_for_review_gifs'
mkdir(expert_review_gif_dir)

cd(expert_review_gif_dir)

for i=1:num_cases_selected
    nn = case_figs_selected{i}
    [p, nn, ext] = fileparts(nn)

    nn

    h = openfig(case_figs_selected{i});
    MR_toolbox
    pause
    close(h)
end

[names, num] = findFILE(expert_review_gif_dir, '*.gif');
A = cell(num, 1);
for i=1:num
    nn = names{i}
    [p, nn, ext] = fileparts(nn)
    A{i, 1} = nn;
end
filename = '/data/raw_data/retro_cine_3T/dicom_final/imaging_transformer_retrocine_snr_level_tests_expert.xlsx';
writecell(A, filename)

for i=1:num
    nn = names{i}
    [p, nn, ext] = fileparts(nn)

    % if exist(fullfile(expert_review_gif_dir, [nn '_ES.fig']))
    %     continue
    % end

    fig_name = fullfile(expert_review_dir, [nn '.fig']);
    h = openfig(fig_name)
    MR_toolbox
        
    disp('Upsample in all figures ...')

    [I, xlim, ylim, clim, position]=getimagedata(h);

    ED = Matlab_gt_resize_2D_image(double(squeeze(I(:,:,2,:))), 4*size(I, 1), 4*size(I, 2), 5);
    ES = Matlab_gt_resize_2D_image(double(squeeze(I(:,:,11,:))), 4*size(I, 1), 4*size(I, 2), 5);

    h1 = figure('Name', [nn '_ED']); imagescn(ED, clim, [2 6], 32);    
    setimagescn(4*xlim, 4*ylim, clim);

    h2 = figure('Name', [nn '_ES']); imagescn(ES, clim, [2 6], 32);
    setimagescn(4*xlim, 4*ylim, clim);

    % h2.WindowState = 'maximized';
    % h1.WindowState = 'maximized';
    
    dx = 4*xlim(1, 1)+15
    dy = 4*mean(ylim(1, 2))-40

    axes = findall(h1,'type','axes');
    text(axes(11), dx, dy, 'GT', 'FontSize', 20, 'Color', [1 1 0], 'FontWeight', 'bold', 'BackgroundColor', [0 0.5 0.5]);
    for k=1:10
        text(axes(k), dx, dy, num2str(10-k+1), 'FontSize', 20, 'Color', [1 1 0], 'FontWeight', 'bold', 'BackgroundColor', [0 0.5 0.5]);
    end

    axes = findall(h2,'type','axes');
    text(axes(11), dx, dy, 'GT', 'FontSize', 20, 'Color', [1 1 0], 'FontWeight', 'bold', 'BackgroundColor', [0 0.5 0.5]);
    for k=1:10
        text(axes(k), dx, dy, num2str(10-k+1), 'FontSize', 20, 'Color', [1 1 0], 'FontWeight', 'bold', 'BackgroundColor', [0 0.5 0.5]);
    end

    disp('save all figures ...')
    %pause

    saveas(h1, fullfile(expert_review_gif_dir, [nn '_ED.bmp']), 'bmp')
    saveas(h2, fullfile(expert_review_gif_dir, [nn '_ES.bmp']), 'bmp')
    % saveas(h1, fullfile(expert_review_dir, [nn '_ED.fig']), 'fig')
    % saveas(h2, fullfile(expert_review_dir, [nn '_ES.fig']), 'fig')

    disp('close all figures ...')

    closeall
end

% make show cases
expert_review_fig_dir = '/data/raw_data/retro_cine_3T/dicom_final/expert_cases_200m_for_review_fig'
mkdir(expert_review_fig_dir)

ind_snr = find(snr_levels==0.2)

for i=1:num_cases_selected
    nn = case_figs_selected{i}
    [p, nn, ext] = fileparts(nn)

    nn

    if ~isempty(strfind(nn, 'ED')) | ~isempty(strfind(nn, 'ES'))
        continue
    end

    ind = find(nn=='_');

    case_name = nn(1:ind(1)-1)
    view_str = nn(ind(1)+1:ind(2)-1)
    slc_num = nn(ind(3)+1:end);
    slc_num = str2num(slc_num)

    case_dir = fullfile(raw_200m_dir, case_name);

    a = readNPY(fullfile(case_dir, [view_str '_input.npy']));
    size(a)

    noise_sigmas = readNPY(fullfile(case_dir, [view_str '_noise_sigmas.npy']));

    if strcmp(view_str, 'sax')==1
        gt = complex(readNPY(fullfile(case_dir, ['gt_' view_str '_real.npy'])), readNPY(fullfile(case_dir, ['gt_' view_str '_imag.npy'])));
        size(gt)
        gt = abs(gt(:,:,:,slc_num+1));

        b = complex(readNPY(fullfile(case_dir, [view_str '_output_real.npy'])), readNPY(fullfile(case_dir, [view_str '_output_imag.npy'])));
        size(b)
        b = abs(b(:,:,:,slc_num+1,:));
    else

        if strcmp(view_str, 'ch2')
            v_str = '2ch'
        else
            v_str = '4ch'
        end

        gt = readNPY(fullfile(case_dir, ['gt_' v_str '.npy']));
        size(gt)

        b = readNPY(fullfile(case_dir, [view_str '_output.npy']));
        size(b)
    end

    a = squeeze(a);
    b = squeeze(b);
    h = figure('Name', ''); imagescn(cat(4, a(:,:,:,ind_snr), b(:,:,:,ind_snr), gt), [0 3.5*mean(gt(:))], [1 3], [12], 3);

    axes = findall(h,'type','axes');

    RO = size(a, 1)
    E1 = size(a, 2)

    dx = E1/2-15
    dy = RO-10

    text(axes(1), dx, dy, 'GT', 'FontSize', 20, 'Color', [1 1 0], 'FontWeight', 'bold', 'BackgroundColor', [0 0.5 0.5]);
    text(axes(2), E1/2-40, dy, 'Model output', 'FontSize', 20, 'Color', [1 1 0], 'FontWeight', 'bold', 'BackgroundColor', [0 0.5 0.5]);
    text(axes(3), dx, dy, 'Input', 'FontSize', 20, 'Color', [1 1 0], 'FontWeight', 'bold', 'BackgroundColor', [0 0.5 0.5]);

    saveas(h, fullfile(expert_review_fig_dir, [nn '.fig']), 'fig');
end

%% EF and other biomarkers

res_dir = res_200m_dir
raw_dir = raw_200m_dir
csv_file = 'results_summary_20250308173403.csv'

res_dir = res_100m_dir
raw_dir = raw_200m_dir
csv_file = 'results_summary_20250308173403.csv'

res_dir = res_28m_dir
raw_dir = raw_200m_dir
csv_file = 'results_summary_20250308173403.csv'


res_dict = dictionary;
error_keys = [];

for i=1:num_dirs
    key = cases_dirs{i}

    raw_case_dir = fullfile(raw_dir, key)
    case_dir = fullfile(res_dir, key)

    EDV = zeros(num_levels, 1);
    ESV = EDV;
    EF = EDV;
    MASS = EDV;
    EDV_detailed = EDV;
    ESV_detailed = EDV;
    EF_detailed = EDV;
    MASS_detailed = EDV;

    v = table(EDV, ESV, EF, MASS, EDV_detailed, ESV_detailed, EF_detailed, MASS_detailed);

    [name_csvs, num_csv] = findFILE(case_dir, '*.csv');

    if num_csv==1
        [fpath, fname, ext] = fileparts(name_csvs{1});
        csv_file = [fname ext];
    else
        csv_time = 0;
        csv_file = [];
        for nc = 1:num_csv
            [fpath, fname, ext] = fileparts(name_csvs{nc});
            ind = find(fname == '_');
            a_csv_time = str2num(fname(ind(end)+1:end));
            if a_csv_time > csv_time
                csv_time = a_csv_time;
                csv_file = [fname ext];
            end
        end
    end

    raw = readtable(fullfile(case_dir, csv_file));

    nn = size(raw.pt_id, 1);
    for n=1:nn
        pt = raw.pt_id{n};
        num_str = pt(5:end);
        if strcmp(num_str, 'GT')
            snr_ind = num_levels+1;
        else
            snr_ind = str2num(pt(5:end)) + 1;
        end
        v.EDV(snr_ind) = raw.LV_volume_params_rounded_EDV(n);
        v.ESV(snr_ind) = raw.LV_volume_params_rounded_ESV(n);
        v.EF(snr_ind) = 100*raw.LV_volume_params_rounded_EF(n);
        v.MASS(snr_ind) = raw.LV_mass_params_LVM_rounded(n);
        v.EDV_detailed(snr_ind) = raw.LV_volume_params_detailed_EDV(n);
        v.ESV_detailed(snr_ind) = raw.LV_volume_params_detailed_ESV(n);
        v.EF_detailed(snr_ind) = 100*raw.LV_volume_params_detailed_EF(n);
        v.MASS_detailed(snr_ind) = raw.LV_mass_params_LVM_detailed(n);
    end

    if v.EF(end-1)>0.01
        res_dict(key) = {v};
    else
        error_keys = [error_keys {key}]
        disp(['Results are not computed for  ' key ])
    end
end

error_keys

save(fullfile(res_dir, 'EF_results.mat'), 'res_dict', 'error_keys')

% ---------------------------------------------------------------------

num_d = numEntries(res_dict)
E = entries(res_dict)

keys = [];
EDV = zeros(num_d, num_levels+1);
ESV = zeros(num_d, num_levels+1);
EF = zeros(num_d, num_levels+1);
MASS = zeros(num_d, num_levels+1);
EDV_detailed = zeros(num_d, num_levels+1);
ESV_detailed = zeros(num_d, num_levels+1);
EF_detailed = zeros(num_d, num_levels+1);
MASS_detailed = zeros(num_d, num_levels+1);

for i=1:num_d
    EDV(i, :) = E.Value{i}.EDV;
    ESV(i, :) = E.Value{i}.ESV;
    EF(i, :) = E.Value{i}.EF;
    MASS(i, :) = E.Value{i}.MASS;
    EDV_detailed(i, :) = E.Value{i}.EDV_detailed;
    ESV_detailed(i, :) = E.Value{i}.ESV_detailed;
    MASS_detailed(i, :) = E.Value{i}.MASS_detailed;
    EF_detailed(i, :) = E.Value{i}.EF_detailed;    
    keys(i, 1) = E.Key(i);
end

CNR = zeros(num_d, num_levels+1);
for i=1:num_d
    k = E.Key(i)
    res_case_dir = fullfile(raw_dir, k);
    if ~exist(fullfile(res_case_dir, 'record.mat'))
        continue
    end
    record = load(fullfile(res_case_dir, 'record.mat'));
    slc_ind = find(record.cnr_gt_slices(1, :)>0);
    CNR(i, 1:end-1) = mean(record.cnr_input_slices(:, slc_ind), 2);
    CNR(i, end) = mean(record.cnr_gt_slices(1, slc_ind), 2);
end

cine_ai_analysis_res = table(keys, EDV, ESV, EF, MASS, EDV_detailed, ESV_detailed, EF_detailed, MASS_detailed, CNR);

success_rate = zeros(1, num_levels);

for k=1:num_levels
    ind = find(EF(:, k) > 0.01);
    success_rate(k) = 100 * (numel(ind)/num_d);
end

clear xlim ylim
h = figure("Name", "cine-analysis rate");
hold on
plot(snr_levels, success_rate, 'b.', 'MarkerSize', 8);
plot(snr_levels, success_rate, '--', 'LineWidth', 2);
hold off
grid on
box on
xlabel('snr')
ylabel('cine-analysis rate (%)')
xticks(snr_levels)
xscale log
ylim([0 105])

disp([snr_levels; success_rate]')
saveas(h, fullfile(res_dir, 'success_rate.fig'), 'fig');

EF = cine_ai_analysis_res.EF;
diff = abs(cine_ai_analysis_res.EF(:,end-2)-cine_ai_analysis_res.EF(:,end));
diff2 = abs(cine_ai_analysis_res.EF(:,end-2)-cine_ai_analysis_res.EF(:,end));
[v, indd] = sort(diff)
case_ind = find(CNR(:,end)>0 & cine_ai_analysis_res.EF(:, end) > 0.01 & diff<2 & diff2<3);
size(case_ind)

closeall
for kk=1:4
    if kk==1
        v = cine_ai_analysis_res.EF(case_ind, :);
        GT_v = v(:, end);
        y_label_str = "GT-Model, EF (%)"
        x_label_str = "(GT+Model)/2, EF (%)"
        h_fig = figure("Name", "EF, Ground-truth vs. Model")
        y_lim = [-12 12]
        x_lim = [0 100]
        suffix = '%'
    end

    if kk==2
        v = cine_ai_analysis_res.EDV(case_ind, :);
        GT_v = v(:, end);
        y_label_str = "GT-Model, EDV (ml)"
        x_label_str = "(GT+Model)/2, EDV (ml)"
        h_fig = figure("Name", "EDV, Ground-truth vs. Model")
        y_lim = [-40 40]
        x_lim = [0 350]
        suffix = 'ml'
    end

    if kk==3
        v = cine_ai_analysis_res.ESV(case_ind, :);
        GT_v = v(:, end);
        y_label_str = "GT-Model, ESV (ml)"
        x_label_str = "(GT+Model)/2, ESV (ml)"
        h_fig = figure("Name", "ESV, Ground-truth vs. Model")
        y_lim = [-40 40]
        x_lim = [0 350]
        suffix = 'ml'
    end

    if kk==4
        v = cine_ai_analysis_res.MASS(case_ind, :);
        GT_v = v(:, end);
        y_label_str = "GT-Model, MASS (g)"
        x_label_str = "(GT+Model)/2, MASS (g)"
        h_fig = figure("Name", "Myocardium mass, Ground-truth vs. Model")
        y_lim = [-40 40]
        x_lim = [0 350]
        suffix = 'g'
    end
    
    fig_name = fullfile(res_dir, get(h_fig, 'Name'))
    h_fig.WindowState = 'maximized';
    
    abs_mean_diff = zeros(numel(snr_levels), 1);
    CRs = zeros(numel(snr_levels), 2);
    
    symbol = '+'
    markersize = 12
    CR_range = 1.65
    
    for k=num_levels:-1:1
    
        snr_ind = num_levels - k + 1;
    
        snr_levels(snr_ind)
    
        a_v = v(:, k);
        ind = find(a_v>0.01);
    
        B = a_v(ind);
        A = GT_v(ind);
    
        subplot(2, 5, num_levels-k+1)
        hold on
        [means,diffs,meanDiff,CR,linFit] = BlandAltman2(A', B', 2, symbol, 0, markersize, CR_range);
        hold off        
        ylabel(y_label_str)
        xlabel(x_label_str)
        box on
        grid on
        ylim(y_lim)
        xlim(x_lim)
        ha = gca;
        ha.YAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
        ha.YAxis.MinorTickValues = y_lim(1):y_lim(2);
    
        [h, p] = ttest(A', B');
        %title(['SNR = ' num2str(snr_levels(17-k+5+1)) ', mean diff = ' num2str(meanDiff, 4) '%, p=' num2str(p, 5)])    
    
        CRs(snr_ind, :) = CR;
        abs_mean_diff(snr_ind) = mean(abs(A-B));
    
        title({['SNR = ' num2str(snr_levels(k)) ', mean diff = ' num2str(meanDiff, 4) suffix]; ['90% CR = ' num2str(CR(1), 4) ' to ' num2str(CR(2), 4) suffix]})
    end
    
    save([fig_name '_CRs.mat'], 'CRs')
    saveas(h_fig, [fig_name '.fig'], 'fig');
    close(h_fig)
end

% plot with CNR

EF_used = cine_ai_analysis_res.EF(case_ind, 1:end-1);
CNR_used = cine_ai_analysis_res.CNR(case_ind, 1:end-1);

snr_levels_used = snr_levels

size(EF_used)
size(CNR_used)

size(snr_levels_used)

FF = CNR_used(:);
for s=1:132:size(E)
    start_s = s;
    if start_s+132<=size(E,1)
        end_s = start_s+132;
    else
        end_s = size(E);
    end

    mean(Y(start_s:end_s))
    std(Y(start_s:end_s))

    median(Y(start_s:end_s))
end

closeall
CNR_v = cine_ai_analysis_res.CNR(case_ind, 1:end-1);
[Y,E] = sort(CNR_v(:));

size(CNR_v)

for kk=1:4
   
    if kk==1
        v = cine_ai_analysis_res.EF(case_ind, 1:end-1);
        GT_v = cine_ai_analysis_res.EF(case_ind, end);
        y_label_str = "GT-Model, EF (%)"
        x_label_str = "(GT+Model)/2, EF (%)"
        h_fig = figure("Name", "CNR-EF, Groud-truth vs. Model")
        y_lim = [-12 12]
        x_lim = [0 100]
        suffix = '%'
    end

    if kk==2
        v = cine_ai_analysis_res.EDV(case_ind, 1:end-1);
        GT_v = cine_ai_analysis_res.EDV(case_ind, end);
        y_label_str = "GT-Model, EDV (ml)"
        x_label_str = "(GT+Model)/2, EDV (ml)"
        h_fig = figure("Name", "CNR-EDV, Groud-truth vs. Model")
        y_lim = [-40 40]
        x_lim = [0 350]
        suffix = 'ml'
    end

    if kk==3
        v = cine_ai_analysis_res.ESV(case_ind, 1:end-1);
        GT_v = cine_ai_analysis_res.ESV(case_ind, end);
        y_label_str = "GT-Model, ESV (ml)"
        x_label_str = "(GT+Model)/2, ESV (ml)"
        h_fig = figure("Name", "CNR-ESV, Groud-truth vs. Model")
        y_lim = [-40 40]
        x_lim = [0 350]
        suffix = 'ml'
    end

    if kk==4
        v = cine_ai_analysis_res.MASS(case_ind, 1:end-1);
        GT_v = cine_ai_analysis_res.MASS(case_ind, end);
        y_label_str = "GT-Model, MASS (g)"
        x_label_str = "(GT+Model)/2, MASS (g)"
        h_fig = figure("Name", "CNR-MASS, Groud-truth vs. Model")
        y_lim = [-40 40]
        x_lim = [0 350]
        suffix = 'g'
    end
    
    fig_name = fullfile(res_dir, get(h_fig, 'Name'))
    h_fig.WindowState = 'maximized';
    
    CRs = [];
    
    symbol = '+'
    markersize = 12
    CR_range = 1.65
    
    num_bins = 10;
    steps = round(size(E, 1)/num_bins);

    bins = 1:steps:size(E,1);

    CRs = [];
    for s=1:numel(bins)
        start_s = bins(s);
        if s < numel(bins)
            end_s = start_s+steps;
        else
            end_s = size(E, 1);
        end
    
        median_cnr = median(Y(start_s:end_s))

        cnr_ind = E(start_s:end_s, 1);
        v2 = v(:);
        a_v = v2(cnr_ind);
        ind = find(a_v>0.01);
    
        B = a_v(ind);

        [I, J] = ind2sub(size(CNR_v), cnr_ind(ind));

        A = GT_v(I);
    
        subplot(2, ceil(num_bins/2), num_bins-s+1)
        hold on
        [means,diffs,meanDiff,CR,linFit] = BlandAltman2(A', B', 2, symbol, 0, markersize, CR_range)        
        hold off        
        ylabel(y_label_str)
        xlabel(x_label_str)
        box on
        grid on
        ylim(y_lim)
        xlim(x_lim)
        ha = gca;
        ha.YAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
        ha.YAxis.MinorTickValues = y_lim(1):y_lim(2);
    
        [h, p] = ttest(A', B');
        %title(['SNR = ' num2str(snr_levels(17-k+5+1)) ', mean diff = ' num2str(meanDiff, 4) '%, p=' num2str(p, 5)])    
    
        CRs = [CRs; CR];
    
        title({['median CNR = ' num2str(median_cnr, 4) ', mean diff = ' num2str(meanDiff, 4) suffix]; ['90% CR = ' num2str(CR(1), 4) ' to ' num2str(CR(2), 4) suffix]})
    end
    
    save([fig_name '_CRs.mat'], 'CRs')
    saveas(h_fig, [fig_name '.fig'], 'fig');
    close(h_fig)
end

% prepare a seg, cnr, snr plot

cd /data/raw_data/retro_cine_3T/dicom_final/final_200m/00050939/
cd /data/raw_data/retro_cine_3T/dicom_final/final_200m/00050938/
t = res_dict('00050938')
t = t{1}

sax_input = readNPY('sax_input');
size(sax_input)

sax_output_c = complex(readNPY('sax_output_real'), readNPY('sax_output_imag'));
size(sax_output_c)

sax_output = abs(sax_output_c);

seg_mask = readNPY('SNR_GT__rounded_masks.npy');
seg_mask = permute(seg_mask, [3 4 2 1]);
size(seg_mask)

gt_sax_c = complex(readNPY('gt_sax_real'), readNPY('gt_sax_imag'));
gt_sax = abs(gt_sax_c);
size(gt_sax)

figure; imagescn(seg_mask, [], [], [], 3);

h_gt = figure('Name', 'gt'); imagescn(gt_sax(:,:,:,end:-1:1), [0 200], [3 4], [16], 3);
h_input = figure('Name', 'input'); imagescn(sax_input(:,:,:,end:-1:1,4), [0 200], [3 4], [16], 3);
h_output = figure('Name', 'output'); imagescn(sax_output(:,:,:,end:-1:1,4), [0 200], [3 4], [16], 3);

h_mask = figure('Name', 'seg'); imagescn(seg_mask, [], [3 4], [16], 3);

[I, xlim, ylim, clim, position]=getimagedata(h_gt);
figure(h_input)
setimagescn(xlim, ylim);
figure(h_output)
setimagescn(xlim, ylim);
figure(h_mask)
setimagescn(xlim, ylim);

saveas(h_gt, '/data/raw_data/retro_cine_3T/dicom_final/gt_sax.fig', 'fig')
saveas(h_input, '/data/raw_data/retro_cine_3T/dicom_final/sax_input.fig', 'fig')
saveas(h_output, '/data/raw_data/retro_cine_3T/dicom_final/sax_output.fig', 'fig')
saveas(h_mask, '/data/raw_data/retro_cine_3T/dicom_final/seg_mask.fig', 'fig')

%% data statistics

% training
fnames = {'train_3D_3T_retro_cine_2018.h5', 'train_3D_3T_retro_cine_2019.h5', 'train_3D_3T_retro_cine_2020.h5'}

fnames = {'test_cine_sig_32.0_3000.h5'}

fnames = {'retro_cine_3T_20241202.h5'}

data_dir = '/fastdata/mri/'

num_series = zeros(numel(fnames), 1)
f_num_images = zeros(numel(fnames), 1)
for k=1:numel(fnames)

    dd = h5info(fullfile(data_dir, fnames{k}));
    N = numel(dd.Groups)
    num_series(k) = N;

    num_images = 0
    for t=1:N
        if mod(t, 1000) == 0
            disp([num2str(t) ' out of ' num2str(N) ', number of images ' num2str(num_images)])
        end
        a = dd.Groups(t).Name;
        im = h5read(fullfile(data_dir, fnames{k}), [a '/image']);
        if ~isreal(im)
            num_images = num_images + size(im.r, 3);
        else
            num_images = num_images + size(im, 3);
        end
    end
    disp(['number of images ' num2str(num_images)])
    f_num_images(k) = num_images;
end

sum(num_series)
sum(f_num_images)

%% find patient age, mean
data_dir = '/data/raw_data/retro_cine_3T/raw/'

[names, num] = FindSubDirs(data_dir)

pt_ids = []
pt_ages = []
pt_gender = []
SLCs = [];

for n=1:num

    case_dir = fullfile(data_dir, names{n})

    [fnames, K] = findFILE(case_dir, 'Retro_Lin_Cine*GLS*.h5'); 
    for k=1:K
        fnames{k}
        try
            dset = ismrmrd.Dataset(fnames{k});
            hdr = ismrmrd.xml.deserialize(dset.readxml);
            dset.close();
        catch
            continue
        end

        SLCs = [SLCs; hdr.encoding.encodingLimits.slice.maximum+1];

        pt_id = hdr.subjectInformation.patientID;

        found_flag = 0;
        for ii=1:numel(pt_ids)
            if strcmp(pt_id, pt_ids{ii})==1
                found_flag = 1;
                break;
            end
        end

        if (found_flag == 0)
            pt_ids = [pt_ids; {pt_id}];
            pt_gender = [pt_gender; hdr.subjectInformation.patientGender];
    
            for ii=1:numel(hdr.userParameters.userParameterDouble)
                if strcmp('PatientAge', hdr.userParameters.userParameterDouble(ii).name) == 1
                    pt_ages = [pt_ages; hdr.userParameters.userParameterDouble(ii).value];
                end
            end
        end
    end
end

%% create noise training figure 1

cd /data/raw_data/retro_cine_3T/raw_selected/00047882/ch2/Retro_Lin_Cine_2DT_LAX_GLS_000000_591071387_591071396_1285_00000000-000000/res_GT_RetroCine_gmap_augmentation/DebugOutput/

acs = readGTPlusExportData('GFactor_aug_acsSrc_n0_s0_slc0_encoding_0_R5');
size(acs)
figure; plotKSpaceSamplingPattern(acs);
figure; imagescn(abs(acs), [], [3 10])

ker = readGTPlusExportData('GFactor_aug_ker_n0_s0_slc0_encoding_0_R5');
size(ker)

kIm = readGTPlusExportData('GFactor_aug_kIm_n0_s0_slc0_encoding_0_R5');
size(kIm)

cmap = readGTPlusExportData('GFactor_aug_coilMap_n0_s0_slc0_encoding_0_R5');
size(cmap)
figure; imagescn(abs(cmap), [], [3 10])

unmixC = readGTPlusExportData('GFactor_aug_unmixC_n0_s0_slc0_encoding_0_R5');
size(unmixC)
figure; imagescn(abs(unmixC), [], [3 10]);

gmap = readGTPlusExportData('GFactor_aug_gFactor_n0_s0_slc0_encoding_0_R5');
size(gmap)
figure; imagescn(abs(gmap));

gmap_all = analyze75read('gfactor_augmented_encoding_0_1');
size(gmap_all)
figure; imagescn(abs(gmap_all), [], [1 7]);

im = readGTPlusExportData('recon_res_encoding_0_1');
size(im)
figure; imagescn(abs(squeeze(im)), [0 160], [], [], 3);

cd /data/raw_data/retro_cine_3T/raw_selected/00047882/ch2/Retro_Lin_Cine_2DT_LAX_GLS_000000_591071387_591071396_1285_00000000-000000/res_GT_RetroCine_gmap_augmentation/DebugOutput/snr_data/slc_0/

% python3 ./projects/mri_imaging/data/create_test_set_mri_snr_level.py --input_dir /data/raw_data/retro_cine_3T/raw_selected/00047882/ch2/Retro_Lin_Cine_2DT_LAX_GLS_000000_591071387_591071396_1285_00000000-000000/res_GT_RetroCine_gmap_augmentation/DebugOutput --output_dir /data/raw_data/retro_cine_3T/raw_selected/00047882/ch2/Retro_Lin_Cine_2DT_LAX_GLS_000000_591071387_591071396_1285_00000000-000000/res_GT_RetroCine_gmap_augmentation/DebugOutput/snr_data

clean = complex(readNPY('clean_real'), readNPY('clean_imag'));
size(clean)

noisy = complex(readNPY('input_real'), readNPY('input_imag'));
size(noisy)

noise_sigmas = readNPY('noise_sigmas')

ind = 14
noise_sigmas(ind)
figure; imagescn(abs(cat(4, clean(:,:,:,ind), noisy(:,:,:,ind), clean(:,:,:,ind)-noisy(:,:,:,ind))), [], [], [], 3);


%% create figures

ch2_dir = '/data/raw_data/retro_cine_3T/raw_selected/00048236/ch2/Retro_Lin_Cine_2DT_LAX_GLS_000000_640763573_640763582_424_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/'
ch4_dir = '/data/raw_data/retro_cine_3T/raw_selected/00048236/ch4/Retro_Lin_Cine_2DT_LAX_GLS_000000_640763573_640763582_423_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/'
sax_dir = '/data/raw_data/retro_cine_3T/raw_selected/00048236/sax/Retro_Lin_Cine_2DT_LAX_GLS_000000_640763573_640763582_426_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/'

cd /data/raw_data/retro_cine_3T/dicom_final/final_200m/00048236/

ch2_noise_sigmas = readNPY('ch2_noise_sigmas.npy')
ch4_noise_sigmas = readNPY('ch4_noise_sigmas.npy')
sax_noise_sigmas = readNPY('sax_noise_sigmas.npy')

ch2_swin_dir = fullfile(ch2_dir, 'res', 'slc_0', 'hrnet_Swin')
ch4_swin_dir = fullfile(ch4_dir, 'res', 'slc_0', 'hrnet_Swin')
sax_swin_dir = fullfile(sax_dir, 'res', 'slc_4', 'hrnet_Swin')

ch2_vit3D_large_dir = fullfile(ch2_dir, 'res', 'slc_0', 'hrnet_ViT3D_large')
ch4_vit3D_large_dir = fullfile(ch4_dir, 'res', 'slc_0', 'hrnet_ViT3D_large')
sax_vit3D_large_dir = fullfile(sax_dir, 'res', 'slc_4', 'hrnet_ViT3D_large')

ch2_cnnt_large_dir = fullfile(ch2_dir, 'res', 'slc_0', 'hrnet_TTTTTT')
ch4_cnnt_large_dir = fullfile(ch4_dir, 'res', 'slc_0', 'hrnet_TTTTTT')
sax_cnnt_large_dir = fullfile(sax_dir, 'res', 'slc_4', 'hrnet_TTTTTT')

ch2_c3_large_dir = fullfile(ch2_dir, 'res', 'slc_0', 'hrnet_C3_large')
ch4_c3_large_dir = fullfile(ch4_dir, 'res', 'slc_0', 'hrnet_C3_large')
sax_c3_large_dir = fullfile(sax_dir, 'res', 'slc_4', 'hrnet_C3_large')

ch2_swin_output = readNPY(fullfile(ch2_swin_dir, 'output'));
ch2_vit3D_large_output = readNPY(fullfile(ch2_vit3D_large_dir, 'output'));
ch2_cnnt_large_output = complex(readNPY(fullfile(ch2_cnnt_large_dir, 'output_real')), readNPY(fullfile(ch2_cnnt_large_dir, 'output_imag')));
ch2_c3_large_output = readNPY(fullfile(ch2_c3_large_dir, 'output'));

ch4_swin_output = readNPY(fullfile(ch4_swin_dir, 'output'));
ch4_vit3D_large_output = readNPY(fullfile(ch4_vit3D_large_dir, 'output'));
ch4_cnnt_large_output = complex(readNPY(fullfile(ch4_cnnt_large_dir, 'output_real')), readNPY(fullfile(ch4_cnnt_large_dir, 'output_imag')));
ch4_c3_large_output = readNPY(fullfile(ch4_c3_large_dir, 'output'));

sax_swin_output = readNPY(fullfile(sax_swin_dir, 'output'));
sax_vit3D_large_output = readNPY(fullfile(sax_vit3D_large_dir, 'output'));
sax_cnnt_large_output = complex(readNPY(fullfile(sax_cnnt_large_dir, 'output_real')), readNPY(fullfile(sax_cnnt_large_dir, 'output_imag')));
sax_c3_large_output = readNPY(fullfile(sax_c3_large_dir, 'output'));

for k=1:numel(ch2_noise_sigmas)-1
    ch2_swin_output(:,:,:,k) = ch2_swin_output(:,:,:,k) * ch2_noise_sigmas(k);
    ch2_vit3D_large_output(:,:,:,k) = ch2_vit3D_large_output(:,:,:,k) * ch2_noise_sigmas(k);
    ch2_cnnt_large_output(:,:,:,k) = ch2_cnnt_large_output(:,:,:,k) * ch2_noise_sigmas(k);
    ch2_c3_large_output(:,:,:,k) = ch2_c3_large_output(:,:,:,k) * ch2_noise_sigmas(k);
end
for k=1:numel(ch4_noise_sigmas)-1
    ch4_swin_output(:,:,:,k) = ch4_swin_output(:,:,:,k) * ch4_noise_sigmas(k);
    ch4_vit3D_large_output(:,:,:,k) = ch4_vit3D_large_output(:,:,:,k) * ch4_noise_sigmas(k);
    ch4_cnnt_large_output(:,:,:,k) = ch4_cnnt_large_output(:,:,:,k) * ch4_noise_sigmas(k);
    ch4_c3_large_output(:,:,:,k) = ch4_c3_large_output(:,:,:,k) * ch4_noise_sigmas(k);
end
for k=1:size(sax_noise_sigmas,1, 1)-1
    sax_swin_output(:,:,:,k) = sax_swin_output(:,:,:,k) * sax_noise_sigmas(k, 5);
    sax_vit3D_large_output(:,:,:,k) = sax_vit3D_large_output(:,:,:,k) * sax_noise_sigmas(k, 5);
    sax_cnnt_large_output(:,:,:,k) = sax_cnnt_large_output(:,:,:,k) * sax_noise_sigmas(k, 5);
    sax_c3_large_output(:,:,:,k) = sax_c3_large_output(:,:,:,k) * sax_noise_sigmas(k, 5);
end


ch2_input = readNPY('ch2_input.npy');
ch4_input = readNPY('ch4_input.npy');
sax_input = readNPY('sax_input.npy');
size(sax_input)

ch2_output = readNPY('ch2_output.npy');
ch4_output = readNPY('ch4_output.npy');
sax_output = readNPY('sax_output.npy');

gt_2ch = readNPY('gt_2ch.npy');
gt_4ch = readNPY('gt_4ch.npy');
gt_sax = readNPY('gt_sax.npy');

size(sax_output)

cd /data/raw_data/retro_cine_3T/dicom_final/final_28m/00048236/
ch2_28m_output = readNPY('ch2_output.npy');
ch4_28m_output = readNPY('ch4_output.npy');
sax_28m_output = readNPY('sax_output.npy');
cd /data/raw_data/retro_cine_3T/dicom_final/final_50m/00048236/
ch2_50m_output = readNPY('ch2_output.npy');
ch4_50m_output = readNPY('ch4_output.npy');
sax_50m_output = readNPY('sax_output.npy');
cd /data/raw_data/retro_cine_3T/dicom_final/final_100m/00048236/
ch2_100m_output = readNPY('ch2_output.npy');
ch4_100m_output = readNPY('ch4_output.npy');
sax_100m_output = readNPY('sax_output.npy');

cd /data/raw_data/retro_cine_3T/raw_selected/00048236/ch2/Retro_Lin_Cine_2DT_LAX_GLS_000000_640763573_640763582_424_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/res/slc_0/hrnet_200m/
ch2_res_200m = load('metrics.mat')
cd /data/raw_data/retro_cine_3T/raw_selected/00048236/ch2/Retro_Lin_Cine_2DT_LAX_GLS_000000_640763573_640763582_424_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/res/slc_0/hrnet_100m/
ch2_res_100m = load('metrics.mat')
cd /data/raw_data/retro_cine_3T/raw_selected/00048236/ch2/Retro_Lin_Cine_2DT_LAX_GLS_000000_640763573_640763582_424_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/res/slc_0/hrnet_50m/
ch2_res_50m = load('metrics.mat')
cd /data/raw_data/retro_cine_3T/raw_selected/00048236/ch2/Retro_Lin_Cine_2DT_LAX_GLS_000000_640763573_640763582_424_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/res/slc_0/hrnet_28m/
ch2_res_28m = load('metrics.mat')

cd /data/raw_data/retro_cine_3T/raw_selected/00048236/ch4/Retro_Lin_Cine_2DT_LAX_GLS_000000_640763573_640763582_423_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/res/slc_0/hrnet_200m/
ch4_res_200m = load('metrics.mat')
cd /data/raw_data/retro_cine_3T/raw_selected/00048236/ch4/Retro_Lin_Cine_2DT_LAX_GLS_000000_640763573_640763582_423_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/res/slc_0/hrnet_100m/
ch4_res_100m = load('metrics.mat')
cd /data/raw_data/retro_cine_3T/raw_selected/00048236/ch4/Retro_Lin_Cine_2DT_LAX_GLS_000000_640763573_640763582_423_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/res/slc_0/hrnet_50m/
ch4_res_50m = load('metrics.mat')
cd /data/raw_data/retro_cine_3T/raw_selected/00048236/ch4/Retro_Lin_Cine_2DT_LAX_GLS_000000_640763573_640763582_423_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/res/slc_0/hrnet_28m/
ch4_res_28m = load('metrics.mat')

cd /data/raw_data/retro_cine_3T/raw_selected/00048236/sax/Retro_Lin_Cine_2DT_LAX_GLS_000000_640763573_640763582_426_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/res/slc_4/hrnet_200m/
sax_res_200m = load('metrics.mat')
cd /data/raw_data/retro_cine_3T/raw_selected/00048236/sax/Retro_Lin_Cine_2DT_LAX_GLS_000000_640763573_640763582_426_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/res/slc_4/hrnet_100m/
sax_res_100m = load('metrics.mat')
cd /data/raw_data/retro_cine_3T/raw_selected/00048236/sax/Retro_Lin_Cine_2DT_LAX_GLS_000000_640763573_640763582_426_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/res/slc_4/hrnet_50m/
sax_res_50m = load('metrics.mat')
cd /data/raw_data/retro_cine_3T/raw_selected/00048236/sax/Retro_Lin_Cine_2DT_LAX_GLS_000000_640763573_640763582_426_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/res/slc_4/hrnet_28m/
sax_res_28m = load('metrics.mat')

cd /data/raw_data/retro_cine_3T/raw_selected/00048236/ch2/Retro_Lin_Cine_2DT_LAX_GLS_000000_640763573_640763582_424_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/
gt_ch2_c= complex(readNPY(fullfile('./res_R2/', 'output_real')), readNPY(fullfile('./res_R2/', 'output_imag')));

ch2_input_c= complex(readNPY(fullfile('./res/slc_0/', 'input_real')), readNPY(fullfile('./res/slc_0/', 'input_imag')));

h = figure;
imagescn(cat(4, abs(ch2_input), abs(ch2_output), abs(ch2_100m_output), abs(ch2_50m_output), abs(ch2_28m_output)), [0 4*mean(gt_2ch(:))], [5 numel(ch2_noise_sigmas)-1], [32], 3);
saveas(h, fullfile('/data/raw_data/retro_cine_3T/dicom_final', 'ch2_IT.fig'))

h = figure;
imagescn(cat(4, abs(ch4_input), abs(ch4_output), abs(ch4_100m_output), abs(ch4_50m_output), abs(ch4_28m_output)), [0 4*mean(gt_4ch(:))], [5 numel(ch4_noise_sigmas)-1], [32], 3);
saveas(h, fullfile('/data/raw_data/retro_cine_3T/dicom_final', 'ch4_IT.fig'))

h = figure;
imagescn(cat(4, abs(sax_input(:,:,:,4,:)), abs(sax_output(:,:,:,4,:)), abs(sax_100m_output(:,:,:,4,:)), abs(sax_50m_output(:,:,:,4,:)), abs(sax_28m_output(:,:,:,4,:))), [0 6*abs(mean(sax_output(:)))], [5 size(sax_noise_sigmas, 1)-1], [32], 3);
saveas(h, fullfile('/data/raw_data/retro_cine_3T/dicom_final', 'sax_IT.fig'))

h = figure;
imagescn(cat(4, abs(ch2_input), abs(ch2_50m_output), abs(ch2_cnnt_large_output), abs(ch2_swin_output), abs(ch2_vit3D_large_output), abs(ch2_c3_large_output)), [0 4*mean(gt_2ch(:))], [6 numel(ch2_noise_sigmas)-1], [32], 3);
saveas(h, fullfile('/data/raw_data/retro_cine_3T/dicom_final', 'ch2_IT_others.fig'))

h = figure;
imagescn(cat(4, abs(ch4_input), abs(ch4_50m_output), abs(ch4_cnnt_large_output), abs(ch4_swin_output), abs(ch4_vit3D_large_output), abs(ch4_c3_large_output)), [0 4*mean(gt_4ch(:))], [6 numel(ch4_noise_sigmas)-1], [32], 3);
saveas(h, fullfile('/data/raw_data/retro_cine_3T/dicom_final', 'ch4_IT_others.fig'))

h = figure;
imagescn(cat(4, squeeze(abs(sax_input(:,:,:,4,:))), squeeze(abs(sax_50m_output(:,:,:,4,:))), abs(sax_cnnt_large_output), abs(sax_swin_output), abs(sax_vit3D_large_output), abs(sax_c3_large_output)), [0 6*mean(sax_output(:))], [6 size(sax_noise_sigmas, 1)-1], [32], 3);
saveas(h, fullfile('/data/raw_data/retro_cine_3T/dicom_final', 'sax_IT_others.fig'))

% 
% for n=1:fig_num
%     [p, fname, ext] = fileparts(fignames{n});
%     disp([num2str(n) ' out of ' num2str(fig_num) ' - ' fname]);
% 
%     for kk=1:4
%         if kk==1
%             final_dir = '/data/raw_data/retro_cine_3T/dicom_final/final_28m/';
%             %mkdir(final_dir)
%             model_str = 'hrnet_28m'
%         end
%         if kk==2
%             final_dir = '/data/raw_data/retro_cine_3T/dicom_final/final_50m/';
%             %mkdir(final_dir)
%             model_str = 'hrnet_50m'
%         end
%         if kk==3
%             final_dir = '/data/raw_data/retro_cine_3T/dicom_final/final_100m/';
%             %mkdir(final_dir)
%             model_str = 'hrnet_100m'
%         end
%         if kk==4
%             final_dir = '/data/raw_data/retro_cine_3T/dicom_final/final_200m/';
%             %mkdir(final_dir)
%             model_str = 'hrnet_200m'
%         end    
% 
%         ind = find(fname=='_');
% 
%         case_name = fname(1:ind(1)-1)
% 
%         final_case_dir = fullfile(final_dir, case_name)
%         mkdir(final_case_dir)
% 
%         if exist(fullfile(final_case_dir, 'sax_output_real.npy'))
%             delete(fullfile(final_case_dir, 'sax_output_real.npy'));
%             delete(fullfile(final_case_dir, 'sax_output_imag.npy'));
%         end
%     end
% end

%% paper 2 figures

snr_levels = [0.05, 0.1, 0.2, 0.5, 0.75, 1.0, 1.5, 2.0, 4.0, 8.0]

cd /data/raw_data/retro_cine_3T/dicom_final/final_200m/00049403/

a = readNPY('gt_sax.npy');
size(a)

cd /data/raw_data/retro_cine_3T/raw_selected/00049403/sax/Retro_Lin_Cine_2DT_LAX_GLS_000000_941869606_941869615_37_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/

gmaps = readNPY('gmap');
size(gmaps)

t = complex(readNPY(fullfile('.', 'input_real')), readNPY(fullfile('.', 'input_imag')));size(a)

cd /data/raw_data/retro_cine_3T/raw_selected/00049403/sax/Retro_Lin_Cine_2DT_LAX_GLS_000000_941869606_941869615_37_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/res/slc_7/

b = complex(readNPY(fullfile('.', 'input_real')), readNPY(fullfile('.', 'input_imag')));size(b)
noise_sigmas = readNPY('noise_sigmas')
N = size(noise_sigmas, 1)

noisy = b;
for n=1:N-1
    noisy(:,:,:,n) = b(:,:,:,n) * noise_sigmas(n);
end

figure; imagescn(cat(4, a(:,:,:,8), abs(noisy)), [], [1 11], [12], 3);

figure; imagescn(gmaps(:,:,8))

gmap = readNPY('gmap');
figure; imagescn(gmap(:,:,1))

roi = load('/data/raw_data/retro_cine_3T/roi.mat')

[J,BW_endo] = roifill(a(:,:,1), roi.ROI_InfoTable(1,1).ROI_x_original, roi.ROI_InfoTable(1,1).ROI_y_original);
[J,BW] = roifill(a(:,:,1), roi.ROI_InfoTable(2,1).ROI_x_original, roi.ROI_InfoTable(2,1).ROI_y_original);

BW_epi = BW - BW_endo;

figure; imagescn(cat(3, a(:,:,1,8), BW_endo, BW_epi))

ind_endo = find(BW_endo(:)>0);
ind_epi = find(BW_epi(:)>0);
G = gmap(:,:,1);
median(G(:))

raw = a(:,:,1,8);

snr_raw = raw ./G;

figure; imagescn(snr_raw)

cnr_raw = mean(snr_raw(ind_endo)) - mean(snr_raw(ind_epi))

signal_bp_raw = mean(raw(ind_endo))
signal_myo_raw = mean(raw(ind_epi))

cnr_im = zeros(10, 1)
signal_bp = zeros(10, 1)
signal_myo = zeros(10, 1)

for n=1:10
    im = b(:,:,1,n);
    snr_im = im ./G;

    cnr_im(n) = abs(mean(snr_im(ind_endo)) - mean(snr_im(ind_epi)));

    signal_bp(n) = abs(mean(im(ind_endo)));
    signal_myo(n) = abs(mean(im(ind_epi)));
end

h = figure('Name', 'Siganl and CNR')

h1 = subplot(2,1,1)
hold on
plot(snr_levels, signal_bp, 'r', 'LineWidth', 3.0)
plot(snr_levels, signal_myo, 'b', 'LineWidth', 3.0)
plot(snr_levels, ones(10,1)*signal_bp_raw, 'r--', 'LineWidth', 3.0)
plot(snr_levels, ones(10,1)*signal_myo_raw, 'b--', 'LineWidth', 3.0)
plot(snr_levels, signal_bp, 'r.', 'MarkerSize', 2.0)
plot(snr_levels, signal_myo, 'b.', 'MarkerSize', 2.0)
hold off
grid on
box on
legend('blood pool', 'myocardium', 'original, blood pool', 'original, myocardium')
xlabel('median snr')
ylabel('Signal levels')
ylim([0 56])
xlim([0 8])
xticks(snr_levels)
xscale log
h1.FontSize = 16

h2 = subplot(2,1,2)
hold on
plot(snr_levels, cnr_im, 'k', 'LineWidth', 3.0)
plot(snr_levels, ones(10,1)*cnr_raw, 'k--', 'LineWidth', 3.0)
plot(snr_levels, cnr_im, 'k.', 'MarkerSize', 2.0)
hold off
grid on
box on
legend('CNR, blood pool - myocardium', 'original CNR')
xlabel('median snr')
ylabel('CNR, blood pool vs. myocardium')
ylim([0 16])
xlim([0 8])
xticks(snr_levels)
xscale log
h2.FontSize = 16

% create more examples

[v, sorted_ind] = sort(case_res.ssim_200m(:, 2))
case_res.c_list(sorted_ind)

for k=1:3000

    c_dir = case_res.c_list(sorted_ind(end-k+1))
    c_dir = c_dir{1};

    ind = find(c_dir=='/');
    slc_str = c_dir(ind(end)+1:end)

    slc = str2num(slc_str(5:end))+1

    view_str =  c_dir(ind(6)+1:ind(7)-1)
    % 
    % if strcmp(view_str, 'sax') == 0
    %     continue
    % end

    input = complex(readNPY(fullfile(c_dir, 'input_real')), readNPY(fullfile(c_dir, 'input_imag')));
    size(input)


    a = complex(readNPY(fullfile(c_dir, 'hrnet_200m', 'output_real')), readNPY(fullfile(c_dir, 'hrnet_200m', 'output_imag')));
    size(a)

    gt = complex(readNPY(fullfile(c_dir, '../../res_R2', 'output_real')), readNPY(fullfile(c_dir, '../../res_R2', 'output_imag')));
    size(gt)

    gt_slc = gt(:,:,:, slc);

    noise_sigmas = readNPY(fullfile(c_dir, 'noise_sigmas'))

    for n=1:10
        input(:,:,:,n) = input(:,:,:,n) * noise_sigmas(n);
        a(:,:,:,n) = a(:,:,:,n) * noise_sigmas(n);
    end

    dd = cat(4, abs(input(:,:,:,[2 3 4 8])), abs(a(:,:,:,[2 3 4 8])), abs(a(:,:,:,[2 3 4 8])));
    size(dd)

    RO = size(dd, 1)
    E1 = size(dd, 2)
    PHS = size(dd, 3)

    res = zeros(RO, E1, PHS, 15);
    res(:,:,:,[2 3 4 5]) = abs(input(:,:,:,[8 4 3 2]));
    res(:,:,:,6) = gt_slc;
    res(:,:,:,[7 8 9 10]) = abs(a(:,:,:,[8 4 3 2]));
    res(:,:,:,11) = gt_slc;
    res(:,:,:,[12 13 14 15]) = abs(a(:,:,:,[8 4 3 2]));

    h = figure; imagescn(res, [0 8*abs(mean(gt(:)))], [3 5], [16], 3);
    fname = ['example_' num2str(k) '_' view_str '_' slc_str '.fig']
    saveas(h, fullfile('/data/raw_data/retro_cine_3T/dicom_final/cases/', fname))
    close(h)
end


[v, sorted_ind] = sort(case_res.psnr_200m(:, 2))
case_res.c_list(sorted_ind)

for k=1:3000

    c_dir = case_res.c_list(sorted_ind(end-k+1))
    c_dir = c_dir{1};

    ind = find(c_dir=='/');
    slc_str = c_dir(ind(end)+1:end)

    slc = str2num(slc_str(5:end))+1

    view_str =  c_dir(ind(6)+1:ind(7)-1)
    % 
    % if strcmp(view_str, 'sax') == 0
    %     continue
    % end

    input = complex(readNPY(fullfile(c_dir, 'input_real')), readNPY(fullfile(c_dir, 'input_imag')));
    size(input)


    a = complex(readNPY(fullfile(c_dir, 'hrnet_200m', 'output_real')), readNPY(fullfile(c_dir, 'hrnet_200m', 'output_imag')));
    size(a)

    gt = complex(readNPY(fullfile(c_dir, '../../', 'input_real')), readNPY(fullfile(c_dir, '../../', 'input_imag')));
    size(gt)

    gt_slc = gt(:,:,:, slc);

    noise_sigmas = readNPY(fullfile(c_dir, 'noise_sigmas'))

    for n=1:10
        input(:,:,:,n) = input(:,:,:,n) * noise_sigmas(n);
        a(:,:,:,n) = a(:,:,:,n) * noise_sigmas(n);
    end

    dd = cat(4, abs(input(:,:,:,[2 3 4 8])), abs(a(:,:,:,[2 3 4 8])), abs(a(:,:,:,[2 3 4 8])));
    size(dd)

    RO = size(dd, 1)
    E1 = size(dd, 2)
    PHS = size(dd, 3)

    res = zeros(RO, E1, PHS, 15);
    res(:,:,:,[2 3 4 5]) = abs(input(:,:,:,[8 4 3 2]));
    res(:,:,:,6) = gt_slc;
    res(:,:,:,[7 8 9 10]) = abs(a(:,:,:,[8 4 3 2]));
    res(:,:,:,11) = gt_slc;
    res(:,:,:,[12 13 14 15]) = abs(a(:,:,:,[8 4 3 2]));

    h = figure; imagescn(res, [0 12*abs(mean(gt(:)))], [3 5], [24], 3);
    fname = ['example_' num2str(k) '_' view_str '_' slc_str '.fig']
    saveas(h, fullfile('/data/raw_data/retro_cine_3T/dicom_final/cases2/', fname))
    close(h)
end

%% view different R data

slc=0

cd /data/raw_data/retro_cine_3T/raw_selected/00047850/ch2/Retro_Lin_Cine_2DT_LAX_GLS_000000_591071167_591071176_821_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/
cd /data/raw_data/retro_cine_3T/raw_selected/00047850/ch4/Retro_Lin_Cine_2DT_LAX_GLS_000000_591071167_591071176_819_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/
cd /data/raw_data/retro_cine_3T/raw_selected/00047850/ch4/Retro_Lin_Cine_2DT_LAX_GLS_000000_591071167_591071176_819_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/
cd /data/raw_data/retro_cine_3T/raw_selected/00047882/ch2/Retro_Lin_Cine_2DT_LAX_GLS_000000_591071387_591071396_1285_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/
cd /data/raw_data/retro_cine_3T/raw_selected/00047881/ch2/Retro_Lin_Cine_2DT_LAX_GLS_000000_591071360_591071369_1220_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/

cd /data/raw_data/retro_cine_3T/raw_selected/00047851/ch2/Retro_Lin_Cine_2DT_LAX_GLS_000000_591071194_591071203_878_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model

cd /data/raw_data/retro_cine_3T/raw_selected/00047852/ch2/Retro_Lin_Cine_2DT_LAX_GLS_000000_591071221_591071230_926_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/

cd /data/raw_data/retro_cine_3T/raw_selected/00047850/sax/Retro_Lin_Cine_2DT_LAX_GLS_000000_591071167_591071176_825_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/
slc = 4

gt = complex(readNPY('input_real'), readNPY('input_imag'));
size(gt)

RO = size(gt, 1)
E1 = size(gt, 2)

noise_sigmas = readNPY(['./res_R_2/slc_' num2str(slc) '/noise_sigmas']);

gmap_R2 = readNPY(['./res_R_2/slc_' num2str(slc) '/gmap']);
gmap_R3 = readNPY(['./res_R_3/slc_' num2str(slc) '/gmap']);
gmap_R4 = readNPY(['./res_R_4/slc_' num2str(slc) '/gmap']);
gmap_R5 = readNPY(['./res_R_5/slc_' num2str(slc) '/gmap']);

gmap_R2 = gmap_R2(:,:,1);
gmap_R3 = gmap_R3(:,:,1);
gmap_R4 = gmap_R4(:,:,1);
gmap_R5 = gmap_R5(:,:,1);

size(gmap_R5)

% median SNR
PHS = size(gt, 3)
snr_R2 = abs(gt) ./ repmat(gmap_R2, [1 1 PHS]);
snr_R3 = abs(gt) ./ repmat(gmap_R3, [1 1 PHS]);
snr_R4 = abs(gt) ./ repmat(gmap_R4, [1 1 PHS]);
snr_R5 = abs(gt) ./ repmat(gmap_R5, [1 1 PHS]);

N = size(noise_sigmas, 1)-1
input_snrs = zeros(4, N);

for n=1:N    
    input_snrs(1, n) = median(snr_R2(:) / noise_sigmas(n));
    input_snrs(2, n) = median(snr_R3(:) / noise_sigmas(n));
    input_snrs(3, n) = median(snr_R4(:) / noise_sigmas(n));
    input_snrs(4, n) = median(snr_R5(:) / noise_sigmas(n));
end

% input
data_R2 = complex(readNPY(['./res_R_2/slc_' num2str(slc) '/input_real']), readNPY(['./res_R_2/slc_' num2str(slc) '/input_imag']));
data_R3 = complex(readNPY(['./res_R_3/slc_' num2str(slc) '/input_real']), readNPY(['./res_R_3/slc_' num2str(slc) '/input_imag']));
data_R4 = complex(readNPY(['./res_R_4/slc_' num2str(slc) '/input_real']), readNPY(['./res_R_4/slc_' num2str(slc) '/input_imag']));
data_R5 = complex(readNPY(['./res_R_5/slc_' num2str(slc) '/input_real']), readNPY(['./res_R_5/slc_' num2str(slc) '/input_imag']));

size(data_R5)

for k=1:N
    data_R2(:,:,:,k) = data_R2(:,:,:,k) * noise_sigmas(k);
    data_R3(:,:,:,k) = data_R3(:,:,:,k) * noise_sigmas(k);
    data_R4(:,:,:,k) = data_R4(:,:,:,k) * noise_sigmas(k);
    data_R5(:,:,:,k) = data_R5(:,:,:,k) * noise_sigmas(k);
end

% output
output_R2 = complex(readNPY(['./res_R_2/slc_' num2str(slc) '/hrnet_200m/output_real']), readNPY(['./res_R_2/slc_' num2str(slc) '/hrnet_200m/output_imag']));
output_R3 = complex(readNPY(['./res_R_3/slc_' num2str(slc) '/hrnet_200m/output_real']), readNPY(['./res_R_3/slc_' num2str(slc) '/hrnet_200m/output_imag']));
output_R4 = complex(readNPY(['./res_R_4/slc_' num2str(slc) '/hrnet_200m/output_real']), readNPY(['./res_R_4/slc_' num2str(slc) '/hrnet_200m/output_imag']));
output_R5 = complex(readNPY(['./res_R_5/slc_' num2str(slc) '/hrnet_200m/output_real']), readNPY(['./res_R_5/slc_' num2str(slc) '/hrnet_200m/output_imag']));

for k=1:N
    output_R2(:,:,:,k) = output_R2(:,:,:,k) * noise_sigmas(k);
    output_R3(:,:,:,k) = output_R3(:,:,:,k) * noise_sigmas(k);
    output_R4(:,:,:,k) = output_R4(:,:,:,k) * noise_sigmas(k);
    output_R5(:,:,:,k) = output_R5(:,:,:,k) * noise_sigmas(k);
end

output_R2 = double(output_R2);
output_R3 = double(output_R3);
output_R4 = double(output_R4);
output_R5 = double(output_R5);

ssims = zeros(4, N);
for n=1:N    
    ssims(1, n) = multissim3(abs(gt), abs(output_R2(:,:,:,n)));
    ssims(2, n) = multissim3(abs(gt), abs(output_R3(:,:,:,n)));
    ssims(3, n) = multissim3(abs(gt), abs(output_R4(:,:,:,n)));
    ssims(4, n) = multissim3(abs(gt), abs(output_R5(:,:,:,n)));
end

max_w = 3*mean(abs(gt(:)));

h_gt = figure; imagescn(abs(gt(:,:,:,slc+1)), [0 max_w], [], [8], 3);

h_gmap = figure; imagescn(abs(cat(3, gmap_R2, gmap_R3, gmap_R4, gmap_R5)), [0 6], [1 4], [32]); colormap('jet');

g2 = median(gmap_R2(:))
g3 = median(gmap_R3(:))
g4 = median(gmap_R4(:))
g5 = median(gmap_R5(:))

axes = findall(h_gmap,'type','axes');

x = E1/2-55
y = RO-10
text(axes(1), x, y, ['R=5, median ' num2str(g5, 5)], 'FontSize', 24, 'Color', [1 1 0], 'FontWeight', 'bold');
text(axes(2), x, y, ['R=4, median ' num2str(g4, 5)], 'FontSize', 24, 'Color', [1 1 0], 'FontWeight', 'bold');
text(axes(3), x, y, ['R=3, median ' num2str(g3, 5)], 'FontSize', 24, 'Color', [1 1 0], 'FontWeight', 'bold');
text(axes(4), x, y, ['R=2, median ' num2str(g2, 5)], 'FontSize', 24, 'Color', [1 1 0], 'FontWeight', 'bold');

input_snrs

h_input = figure; imagescn(abs(cat(4, data_R2, data_R3, data_R4, data_R5)), [0 max_w], [4 10], [32], 3);
axes = findall(h_input,'type','axes');

text(axes(10), 0, y, 'R=5', 'FontSize', 20, 'Color', [1 1 0], 'FontWeight', 'bold', 'BackgroundColor', [0 0.5 0.5]);
text(axes(20), 0, y, 'R=4', 'FontSize', 20, 'Color', [1 1 0], 'FontWeight', 'bold', 'BackgroundColor', [0 0.5 0.5]);
text(axes(30), 0, y, 'R=3', 'FontSize', 20, 'Color', [1 1 0], 'FontWeight', 'bold', 'BackgroundColor', [0 0.5 0.5]);
text(axes(40), 0, y, 'R=2', 'FontSize', 20, 'Color', [1 1 0], 'FontWeight', 'bold', 'BackgroundColor', [0 0.5 0.5]);

for k=1:10
    text(axes(k), x, y, ['median snr ' num2str(input_snrs(4, 10-k+1), 3)], 'FontSize', 16, 'Color', [1 1 1], 'FontWeight', 'bold', 'BackgroundColor', [0 0 0]);
    text(axes(10+k), x, y, ['median snr ' num2str(input_snrs(3, 10-k+1), 3)], 'FontSize', 16, 'Color', [1 1 1], 'FontWeight', 'bold', 'BackgroundColor', [0 0 0]);
    text(axes(20+k), x, y, ['median snr ' num2str(input_snrs(2, 10-k+1), 3)], 'FontSize', 16, 'Color', [1 1 1], 'FontWeight', 'bold', 'BackgroundColor', [0 0 0]);
    text(axes(30+k), x, y, ['median snr ' num2str(input_snrs(1, 10-k+1), 3)], 'FontSize', 16, 'Color', [1 1 1], 'FontWeight', 'bold', 'BackgroundColor', [0 0 0]);
end

h_output = figure; imagescn(abs(cat(4, output_R2, output_R3, output_R4, output_R5)), [0 max_w], [4 10], [32], 3);
axes = findall(h_output,'type','axes');
text(axes(10), 0, y, 'R=5', 'FontSize', 20, 'Color', [1 1 0], 'FontWeight', 'bold', 'BackgroundColor', [0 0.5 0.5]);
text(axes(20), 0, y, 'R=4', 'FontSize', 20, 'Color', [1 1 0], 'FontWeight', 'bold', 'BackgroundColor', [0 0.5 0.5]);
text(axes(30), 0, y, 'R=3', 'FontSize', 20, 'Color', [1 1 0], 'FontWeight', 'bold', 'BackgroundColor', [0 0.5 0.5]);
text(axes(40), 0, y, 'R=2', 'FontSize', 20, 'Color', [1 1 0], 'FontWeight', 'bold', 'BackgroundColor', [0 0.5 0.5]);


for k=1:10
    text(axes(k), x, y, ['SSIM ' num2str(ssims(4, 10-k+1), 3)], 'FontSize', 16, 'Color', [1 1 1], 'FontWeight', 'bold', 'BackgroundColor', [0 0 0]);
    text(axes(10+k), x, y, ['SSIM ' num2str(ssims(3, 10-k+1), 3)], 'FontSize', 16, 'Color', [1 1 1], 'FontWeight', 'bold', 'BackgroundColor', [0 0 0]);
    text(axes(20+k), x, y, ['SSIM ' num2str(ssims(2, 10-k+1), 3)], 'FontSize', 16, 'Color', [1 1 1], 'FontWeight', 'bold', 'BackgroundColor', [0 0 0]);
    text(axes(30+k), x, y, ['SSIM ' num2str(ssims(1, 10-k+1), 3)], 'FontSize', 16, 'Color', [1 1 1], 'FontWeight', 'bold', 'BackgroundColor', [0 0 0]);
end

saveas(h_gmap, fullfile('/data/raw_data/retro_cine_3T/dicom_final/multi_R_gmap.fig'), 'fig');
saveas(h_input, fullfile('/data/raw_data/retro_cine_3T/dicom_final/multi_R_input.fig'), 'fig');
saveas(h_output, fullfile('/data/raw_data/retro_cine_3T/dicom_final/multi_R_output.fig'), 'fig');
saveas(h_gt, fullfile('/data/raw_data/retro_cine_3T/dicom_final/multi_R_gt.fig'), 'fig');

%% first show case

% 00049985_sax_slc_5


cd /data/raw_data/retro_cine_3T/raw_selected/00049985/sax/Retro_Lin_Cine_2DT_LAX_GLS_000000_006537362_006537371_232_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/
slc = 6

gt = complex(readNPY('input_real'), readNPY('input_imag'));
size(gt)

RO = size(gt, 1)
E1 = size(gt, 2)

figure; imagescn(abs(gt(:,:,:,slc)), [], [], [], 3)

noise_sigmas = readNPY('./res/slc_5/noise_sigmas');
gmap_R4 = readNPY(['./res/slc_5/gmap']);

noisy = complex(readNPY('./res/slc_5/input_real'), readNPY('./res/slc_5/input_imag'));
res = complex(readNPY('./res/slc_5/hrnet_200m/output_real'), readNPY('./res/slc_5/hrnet_200m/output_imag'));

for i=1:10
    noisy(:,:,:,i) = noisy(:,:,:,i) * noise_sigmas(i);
    res(:,:,:,i) = res(:,:,:,i) * noise_sigmas(i);
end

ind = 10 -7 + 1
snr_levels(ind)

PHS = size(gt, 3)

closeall
h = figure('Name', '00049985_sax_slc_5')
A = abs(gt(:,:,:,slc));
B = repmat(gmap_R4(:,:,1), [1 1 PHS]);
C = abs(noisy(:,:,:,ind));
D = abs(res(:,:,:,ind));
imagescn(cat(4, C, B, D, A), [], [1 4], [32], 3);

%% compute power spectra
cd /data/raw_data/retro_cine_3T/power_spectra/

RO = 256
E1 = 192

f_pf1_RO = Matlab_gt_generate_kspace_filter_asymmetric(RO, RO*2/8, RO-1, 12, 0);
f_pf2_RO = Matlab_gt_generate_kspace_filter_asymmetric(RO, RO*2/8, RO-1, 24, 0);
f_pf3_RO = Matlab_gt_generate_kspace_filter_asymmetric(RO, RO*2/8, RO-1, 48, 0);

f_pf1_E1 = Matlab_gt_generate_kspace_filter_asymmetric(E1, E1*2/8, E1-1, 12, 0);
f_pf2_E1 = Matlab_gt_generate_kspace_filter_asymmetric(E1, E1*2/8, E1-1, 24, 0);
f_pf3_E1 = Matlab_gt_generate_kspace_filter_asymmetric(E1, E1*2/8, E1-1, 48, 0);

f1_RO = Matlab_gt_generate_kspace_filter_symmetric(RO, 'Gaussian', 1.0, 0);
f2_RO = Matlab_gt_generate_kspace_filter_symmetric(RO, 'Gaussian', 1.5, 0);
f3_RO = Matlab_gt_generate_kspace_filter_symmetric(RO, 'Gaussian', 2.5, 0);

f1_E1 = Matlab_gt_generate_kspace_filter_symmetric(E1, 'Gaussian', 1.0, 0);
f2_E1 = Matlab_gt_generate_kspace_filter_symmetric(E1, 'Gaussian', 1.5, 0);
f3_E1 = Matlab_gt_generate_kspace_filter_symmetric(E1, 'Gaussian', 2.5, 0);

f1_RO(1) = f1_RO(2);
f2_RO(1) = f2_RO(2);
f3_RO(1) = f3_RO(2);

f1_E1(1) = f1_E1(2);
f2_E1(1) = f2_E1(2);
f3_E1(1) = f3_E1(2);

% load a gmap
res_dir = '/data/raw_data/retro_cine_3T/raw_selected/00050644/sax/Retro_Lin_Cine_2DT_LAX_GLS_000000_082295846_082295855_45_00000000-000000/res_GT_RetroCine_gmap_augmentation_SCC_for_AI_OFFLINE2/model/'
gmap = readNPY(fullfile(res_dir, 'gmaps_aug'));
gmap = gmap(:,:,4,2);
figure; imagescn(gmap)

figure; hold on; 
plot(f_pf1_E1/max(f_pf1_E1), 'r', 'LineWidth', 2.0); 
plot(f_pf2_E1/max(f_pf2_E1), 'g', 'LineWidth', 2.0); 
plot(f_pf3_E1/max(f_pf3_E1), 'k', 'LineWidth', 2.0);
hold off
box on
grid on
legend('PF 6/8, T12', 'PF 6/8, T24', 'PF 6/8, T48','Location','SouthEast');
ylim([0 1.1])
xlim([0 E1])
h = gca;
h.FontSize = 16

figure; hold on; 
plot(f1_RO/max(f1_RO), 'r', 'LineWidth', 2.0); 
plot(f2_RO/max(f2_RO), 'g', 'LineWidth', 2.0); 
plot(f3_RO/max(f3_RO), 'k', 'LineWidth', 2.0);
hold off
box on
grid on
legend('Gaussian 1.0', 'Gaussian 1.5', 'Gaussian 2.5','Location','NorthEast');
ylim([0 1.1])
xlim([0 RO])
h = gca;
h.FontSize = 16

% test filters

REP = 256;

nn = complex(randn(RO, E1, REP), randn(RO, E1, REP));
nn_g = nn .* repmat(gmap, [1 1 REP]);

kspace= fft2c(nn_g);

d3 = performKSpaceFilter2D(kspace, f3_RO, f3_E1);
d3 = ifft2c(d3);

figure; imagescn(abs(cat(4, nn, nn_g, d3)), [], [], [], 3);

d4 = performKSpaceFilter2D(kspace, f3_RO, f_pf3_E1);
d4 = ifft2c(d4);

h = figure; imagescn(abs(cat(4, nn, nn_g, d3, d4)), [], [1 4], [], 3);
saveas(h, 'noise_pattern.fig', 'fig')

% compute power spectrum and autocorrelation
d1 = performKSpaceFilter2D(kspace, f1_RO, f1_E1);
d2 = performKSpaceFilter2D(kspace, f2_RO, f2_E1);
d3 = performKSpaceFilter2D(kspace, f3_RO, f3_E1);
d4 = performKSpaceFilter2D(kspace, f3_RO, f_pf3_E1);

fdn = mean(abs(fft2c(nn)), 3);
fdg = mean(abs(fft2c(nn_g)), 3);
fdd1 = mean(abs(d1), 3);
fdd2 = mean(abs(d2), 3);
fdd3 = mean(abs(d3), 3);
fdd4 = mean(abs(d4), 3);

K = RO/2
dn = fdn(K, :); dn = dn / mean(dn(E1/2-1:E1/2+1));
dg = fdg(K, :); dg = dg / mean(dg(E1/2-1:E1/2+1));
dd1 = fdd1(K, :); dd1 = dd1 / mean(dd1(E1/2-1:E1/2+1));
dd2 = fdd2(K, :); dd2 = dd2 / mean(dd2(E1/2-1:E1/2+1));
dd3 = fdd3(K, :); dd3 = dd3 / mean(dd3(E1/2-1:E1/2+1));
dd4 = fdd4(K, :); dd4 = dd4 / mean(dd4(E1/2-1:E1/2+1));

h1 = figure; hold on; 
plot(dn, 'c', 'LineWidth', 4.0); 
plot(dg, 'b', 'LineWidth', 2.0); 
plot(dd1, 'r', 'LineWidth', 2.0); 
plot(dd2, 'g', 'LineWidth', 2.0); 
plot(dd3, 'k', 'LineWidth', 2.0);
plot(dd4, 'm', 'LineWidth', 2.0);
hold off
box on
grid on
legend('White noise', 'with g-factor map', 'Gaussian 1.0', 'Gaussian 1.5', 'Gaussian 2.5','PF6/8 T48','Location','NorthEast');
ylim([0 1.1])
xlim([0 E1])
h = gca;
h.FontSize = 16
saveas(h1, 'power_spectrum.fig', 'fig')

h2 = figure; hold on; 
plot(abs(fftc(dn(:),1))/max((abs(fftc(dn(:),1)))), 'c', 'LineWidth', 4.0); 
plot(abs(fftc(dg(:),1))/max((abs(fftc(dg(:),1)))), 'b', 'LineWidth', 2.0); 
plot(abs(fftc(dd1(:),1))/max((abs(fftc(dd1(:),1)))), 'r', 'LineWidth', 2.0); 
plot(abs(fftc(dd2(:),1))/max((abs(fftc(dd2(:),1)))), 'g', 'LineWidth', 2.0); 
plot(abs(fftc(dd3(:),1))/max((abs(fftc(dd3(:),1)))), 'k', 'LineWidth', 2.0);
plot(abs(fftc(dd4(:),1))/max((abs(fftc(dd4(:),1)))), 'm', 'LineWidth', 2.0);
hold off
box on
grid on
legend('White noise', 'with g-factor map', 'Gaussian 1.0', 'Gaussian 1.5', 'Gaussian 2.5','PF6/8 T48','Location','NorthEast');
xlim([E1/2-8 E1/2+8])
h = gca;
h.FontSize = 16
saveas(h2, 'auto_correlation.fig', 'fig')

%% temp
for k=1:num_cases

    case_dir = fullfile(data_dir, names_cases{k});

   
    [names, num] = FindSubDirs(case_dir);
    
    for i=1:num
   
        if any([(strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0) (strfind(names{i}, 'noise_') > 0)])
            fullfile(case_dir, names{i})
            rmdir(fullfile(case_dir, names{i}), 's')
        end
    end
end

%% db lge

data_dir = '/fastdata/db_lge_data2/raw/'
[names_cases, num_cases] = FindSubDirs(data_dir)

data_dir = '/data/raw_data/db_lge_data/raw/'
[names_cases, num_cases] = FindSubDirs(data_dir)

mat_dir = '/fastdata/db_lge_data2/res/';
moco_key_frame = 4

for k=1:num_cases

    case_dir = fullfile(data_dir, names_cases{k});

    findAndMoveMeasDat(case_dir)
    
    UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;
    
    [names, num] = FindSubDirs(case_dir);
    
    for i=1:num
        UTCases{1, 1} = case_dir;
        UTCases{1, 2} = names{i};
    
        dst_dir = fullfile(mat_dir, names_cases{k});
        if exist(fullfile(dst_dir, [names{i} '_PSIR.mat']))
            continue
        end

        for p=4:5
            if p==1
                UTCases{1, 4} = 'GTPrep_DB_LGE_STCNNT.xml';
                UTCases{1, 5} = 'res_GTPrep_DB_LGE_STCNNT';
            end
            if p==2
                UTCases{1, 4} = 'GTPrep_DB_LGE.xml';
                UTCases{1, 5} = 'res_GTPrep_DB_LGE'; 
            end
            if p==3
                UTCases{1, 4} = 'GTPrep_DB_LGE_half.xml';
                UTCases{1, 5} = 'res_GTPrep_DB_LGE_half'; 
            end
            if p==4
                UTCases{1, 4} = 'GTPrep_DB_LGE_STCNNT_half.xml';
                UTCases{1, 5} = 'res_GTPrep_DB_LGE_STCNNT_half'; 
            end
            if p==5
                UTCases{1, 4} = 'GTPrep_DB_LGE_STCNNT_quarter.xml';
                UTCases{1, 5} = 'res_GTPrep_DB_LGE_STCNNT_quarter'; 
            end
            if p==6
                UTCases{1, 4} = 'GTPrep_DB_LGE_quarter.xml';
                UTCases{1, 5} = 'res_GTPrep_DB_LGE_quarter'; 
            end
            
            fnames = findFILE(fullfile(case_dir, names{i}, UTCases{1, 5}), '*.hdr');
            if (numel(fnames)>0)
                continue
            end       
            
            try
                performUTValidation(UTCases, 0, 0, 'localhost', '9004', 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput');
        
                % debug_dir = fullfile(case_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
                % [data, gmap, res] = prepare_for_stcnnt_inference_LGE(debug_dir); 
            catch
            end
            % if ~isempty(res)
            %     SLC = size(data, 5)
            %     if SLC==1
            %         h = figure('Name', unique_h5names{c}); imagescn(cat(4, data, res), [], [SLC 4], [12], 3);
            %     else
            %         h = figure('Name',unique_h5names{c}); imagescn(cat(4, data, res), [], [4 SLC], [12], 3);
            %     end
            %     % saveas(h, fullfile(data_dir, [names{i} '_' UTCases{1, 5} '.fig']));
            %     close(h)
            % end
        end

        disp([num2str(k) ' - ' case_dir ' - ' names{i}])
    
        % if exist(fullfile(mat_dir, case_dirs{c, 1}, [h5names{c, 1} '--' UTCases{1, 5} '_PSIR.mat']))
        %     continue
        % end
    
        im = [];
        
    
        try
            
    
            [psir_ai_no_ave, header_psir_ai_no_ave, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(case_dir, names{i}, 'res_GTPrep_DB_LGE_STCNNT'), 107, 1);
            [psir_ai_half_no_ave, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(case_dir, names{i}, 'res_GTPrep_DB_LGE_STCNNT_half'), 107, 1);
            [psir_ai_quarter_no_ave, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(case_dir, names{i}, 'res_GTPrep_DB_LGE_STCNNT_quarter'), 107, 1);
        
            [psir_no_ave, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(case_dir, names{i}, 'res_GTPrep_DB_LGE'), 107, 1);
            [psir_half_no_ave, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(case_dir, names{i}, 'res_GTPrep_DB_LGE_half'), 107, 1);
            [psir_quarter_no_ave, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(case_dir, names{i}, 'res_GTPrep_DB_LGE_quarter'), 107, 1);
        
            size(psir_ai_no_ave)
        
            [psir_ai, header_psir_ai, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(case_dir, names{i}, 'res_GTPrep_DB_LGE_STCNNT'), 111, 1);
            [psir_ai_half, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(case_dir, names{i}, 'res_GTPrep_DB_LGE_STCNNT_half'), 111, 1);
            [psir_ai_quarter, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(case_dir, names{i}, 'res_GTPrep_DB_LGE_STCNNT_quarter'), 111, 1);
        
            [psir, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(case_dir, names{i}, 'res_GTPrep_DB_LGE'), 111, 1);
            [psir_half, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(case_dir, names{i}, 'res_GTPrep_DB_LGE_half'), 111, 1);
            [psir_quarter, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(case_dir, names{i}, 'res_GTPrep_DB_LGE_quarter'), 111, 1);
        
            size(psir_ai)
        
            SLC = size(psir_ai, 3);
        
            % debug_dir = fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_STCNNT', 'DebugOutput');
            % [data, gmap, res] = prepare_for_stcnnt_inference_LGE(debug_dir); 
            % 
            % size(data)
            % 
            % raw_psir = squeeze(data(:,:,:,1,:) .* conj(data(:,:,:,2,:)./abs(data(:,:,:,2,:))));
            % raw_psir_ai = squeeze(res(:,:,:,1,:) .* conj(res(:,:,:,2,:)./abs(res(:,:,:,2,:))));
            
            mkdir(dst_dir);
        
            scales = [];
            for slc=1:SLC
                min_psir_ai = header_psir_ai(slc).window_center-header_psir_ai(slc).window_width/2;
                max_psir_ai = header_psir_ai(slc).window_center+header_psir_ai(slc).window_width/2;
                scales = [scales; min_psir_ai max_psir_ai];
            end
    
            h1 = figure('name', [names{i} '-slc' num2str(slc-1)]);
            imagescn(permute(psir, [2 1 3]), [scales], [2 ceil(SLC/2)], [12]);
    
            h_psir_ai = figure('name', [names{i} '-slc' num2str(slc-1)]);
            imagescn(permute(psir_ai, [2 1 3]), [scales], [2 ceil(SLC/2)], [12]);
    
            
            for slc=1:SLC           
            
                im(:,:,:,slc) = cat(3, psir(:,:,slc), psir_ai(:,:,slc), psir_ai_no_ave(:,:,slc, moco_key_frame), ... 
                     psir_half(:,:,slc), psir_ai_half(:,:,slc), psir_ai_half_no_ave(:,:,slc, moco_key_frame), ...
                     psir_quarter(:,:,slc), psir_ai_quarter(:,:,slc), psir_ai_quarter_no_ave(:,:,slc,min([moco_key_frame, size(psir_ai_quarter_no_ave, 4)])) ...
                    );
        
                min_psir_ai = header_psir_ai(slc).window_center-header_psir_ai(slc).window_width/2;
                max_psir_ai = header_psir_ai(slc).window_center+header_psir_ai(slc).window_width/2;
                min_psir_ai_no_ave = header_psir_ai_no_ave(slc).window_center-header_psir_ai_no_ave(slc).window_width/2;
                max_psir_ai_no_ave = header_psir_ai_no_ave(slc).window_center+header_psir_ai_no_ave(slc).window_width/2;
            
                scales = [min_psir_ai max_psir_ai; min_psir_ai max_psir_ai; min_psir_ai max_psir_ai;  min_psir_ai max_psir_ai; min_psir_ai max_psir_ai; min_psir_ai max_psir_ai; ...
                    min_psir_ai_no_ave max_psir_ai_no_ave;  min_psir_ai_no_ave max_psir_ai_no_ave;  min_psir_ai_no_ave max_psir_ai_no_ave];
            
            
                h = figure('name', [names{i} '-slc' num2str(slc-1)]);
                imagescn(permute(im(:,:,:,slc), [2 1 3]), [scales], [3 3], [12]);
                
                %h2 = figure('Name', [h5names{c, 1} '-slc' num2str(slc-1)]); imagescn(real(cat(4, raw_psir(:,:,:,slc), raw_psir_ai(:,:,:,slc))), [-10 10], [1 2], [16], 3);
        
                saveas(h, fullfile(dst_dir, [names{i} '_slc' num2str(slc-1) '_PSIR.fig']), 'fig');
                %saveas(h2, fullfile(dst_dir, [h5names{c, 1} '--' UTCases{1, 5} '_slc' num2str(slc-1) '_PSIR_no_moco_ave.fig']), 'fig');        
            
                close(h)
                %close(h2)
            end  
        
            size(im)
    
            saveas(h1, fullfile(dst_dir, [names{i} '_PSIR.fig']), 'fig');
            saveas(h_psir_ai, fullfile(dst_dir, [names{i} '_PSIR_AI.fig']), 'fig');
            close(h1)
            close(h_psir_ai)
    
            %save(fullfile(dst_dir, [h5names{c, 1} '--' UTCases{1, 5} '_PSIR.mat']), 'im', 'raw_psir', 'raw_psir_ai');
            save(fullfile(dst_dir, [names{i} '_PSIR.mat']), 'im');
        catch
            if exist(dst_dir)
                if size(im, 4) > 2
                    save(fullfile(dst_dir, [names{i} '_PSIR.mat']), 'im');
                else
                    rmdir(dst_dir, 's');
                end
            end
            closeall
        end
    end
end

%% make the mat files
case_dirs = []
h5names = []
ind = 1;
for k=1:num_cases

    case_dir = fullfile(data_dir, names_cases{k});

    findAndMoveMeasDat(case_dir);
    
    UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;
    
    [names, num] = FindSubDirs(case_dir);
    
    for i=1:num
        case_dirs{ind, 1} =  names_cases{k};
        h5names{ind, 1} = names{i};
        ind = ind + 1;
    end
end

numel(h5names)

mat_dir = '/fastdata/db_lge_data2/res/'
mkdir(mat_dir)

moco_key_frame = 3 + 1;

failed_cases = [];

for c=1:numel(case_dirs)

    disp([num2str(c) ' - ' fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1})])

    % if exist(fullfile(mat_dir, case_dirs{c, 1}, [h5names{c, 1} '--' UTCases{1, 5} '_PSIR.mat']))
    %     continue
    % end

    im = [];
    dst_dir = fullfile(mat_dir, case_dirs{c, 1});

    try
        

        [psir_ai_no_ave, header_psir_ai_no_ave, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_STCNNT'), 107, 1);
        [psir_ai_half_no_ave, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_STCNNT_half'), 107, 1);
        [psir_ai_quarter_no_ave, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_STCNNT_quarter'), 107, 1);
    
        [psir_no_ave, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE'), 107, 1);
        [psir_half_no_ave, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_half'), 107, 1);
        [psir_quarter_no_ave, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_quarter'), 107, 1);
    
        size(psir_ai_no_ave)
    
        [psir_ai, header_psir_ai, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_STCNNT'), 111, 1);
        [psir_ai_half, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_STCNNT_half'), 111, 1);
        [psir_ai_quarter, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_STCNNT_quarter'), 111, 1);
    
        [psir, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE'), 111, 1);
        [psir_half, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_half'), 111, 1);
        [psir_quarter, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_quarter'), 111, 1);
    
        size(psir_ai)
    
        SLC = size(psir_ai, 3);
    
        % debug_dir = fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_STCNNT', 'DebugOutput');
        % [data, gmap, res] = prepare_for_stcnnt_inference_LGE(debug_dir); 
        % 
        % size(data)
        % 
        % raw_psir = squeeze(data(:,:,:,1,:) .* conj(data(:,:,:,2,:)./abs(data(:,:,:,2,:))));
        % raw_psir_ai = squeeze(res(:,:,:,1,:) .* conj(res(:,:,:,2,:)./abs(res(:,:,:,2,:))));
        
        mkdir(dst_dir);
    
        scales = [];
        for slc=1:SLC
            min_psir_ai = header_psir_ai(slc).window_center-header_psir_ai(slc).window_width/2;
            max_psir_ai = header_psir_ai(slc).window_center+header_psir_ai(slc).window_width/2;
            scales = [scales; min_psir_ai max_psir_ai];
        end

        h1 = figure('name', [h5names{c, 1} '-slc' num2str(slc-1)]);
        imagescn(permute(psir, [2 1 3]), [scales], [2 ceil(SLC/2)], [12]);

        h_psir_ai = figure('name', [h5names{c, 1} '-slc' num2str(slc-1)]);
        imagescn(permute(psir_ai, [2 1 3]), [scales], [2 ceil(SLC/2)], [12]);

        
        for slc=1:SLC           
        
            im(:,:,:,slc) = cat(3, psir(:,:,slc), psir_ai(:,:,slc), psir_ai_no_ave(:,:,slc, moco_key_frame), ... 
                 psir_half(:,:,slc), psir_ai_half(:,:,slc), psir_ai_half_no_ave(:,:,slc, moco_key_frame), ...
                 psir_quarter(:,:,slc), psir_ai_quarter(:,:,slc), psir_ai_quarter_no_ave(:,:,slc,min([moco_key_frame, size(psir_ai_quarter_no_ave, 4)])) ...
                );
    
            min_psir_ai = header_psir_ai(slc).window_center-header_psir_ai(slc).window_width/2;
            max_psir_ai = header_psir_ai(slc).window_center+header_psir_ai(slc).window_width/2;
            min_psir_ai_no_ave = header_psir_ai_no_ave(slc).window_center-header_psir_ai_no_ave(slc).window_width/2;
            max_psir_ai_no_ave = header_psir_ai_no_ave(slc).window_center+header_psir_ai_no_ave(slc).window_width/2;
        
            scales = [min_psir_ai max_psir_ai; min_psir_ai max_psir_ai; min_psir_ai max_psir_ai;  min_psir_ai max_psir_ai; min_psir_ai max_psir_ai; min_psir_ai max_psir_ai; ...
                min_psir_ai_no_ave max_psir_ai_no_ave;  min_psir_ai_no_ave max_psir_ai_no_ave;  min_psir_ai_no_ave max_psir_ai_no_ave];
        
        
            h = figure('name', [h5names{c, 1} '-slc' num2str(slc-1)]);
            imagescn(permute(im(:,:,:,slc), [2 1 3]), [scales], [3 3], [12]);
            
            %h2 = figure('Name', [h5names{c, 1} '-slc' num2str(slc-1)]); imagescn(real(cat(4, raw_psir(:,:,:,slc), raw_psir_ai(:,:,:,slc))), [-10 10], [1 2], [16], 3);
    
            saveas(h, fullfile(dst_dir, [h5names{c, 1} '--' UTCases{1, 5} '_slc' num2str(slc-1) '_PSIR.fig']), 'fig');
            %saveas(h2, fullfile(dst_dir, [h5names{c, 1} '--' UTCases{1, 5} '_slc' num2str(slc-1) '_PSIR_no_moco_ave.fig']), 'fig');        
        
            close(h)
            %close(h2)
        end  
    
        size(im)

        saveas(h1, fullfile(dst_dir, [h5names{c, 1} '--' UTCases{1, 5} '_PSIR.fig']), 'fig');
        saveas(h_psir_ai, fullfile(dst_dir, [h5names{c, 1} '--' UTCases{1, 5} '_PSIR_AI.fig']), 'fig');
        close(h1)
        close(h_psir_ai)

        %save(fullfile(dst_dir, [h5names{c, 1} '--' UTCases{1, 5} '_PSIR.mat']), 'im', 'raw_psir', 'raw_psir_ai');
        save(fullfile(dst_dir, [h5names{c, 1} '--' UTCases{1, 5} '_PSIR.mat']), 'im');
    catch
        if exist(dst_dir)
            if size(im, 4) > 2
                save(fullfile(dst_dir, [h5names{c, 1} '--' UTCases{1, 5} '_PSIR.mat']), 'im');
            else
                rmdir(dst_dir, 's');
            end
        end
        closeall
        failed_cases = [failed_cases; c case_dirs(c, 1)];
    end

end




%% old code



for k=1:num_cases

    case_dir = fullfile(data_dir, names_cases{k});

    findAndMoveMeasDat(case_dir)
    
    UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;
    
    [names, num] = FindSubDirs(case_dir);

    UTCases{1, 4} = 'GTPrep_DB_LGE.xml';
    UTCases{1, 5} = 'res_GTPrep_DB_LGE'; 

    for i=1:num
        if any([(strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0) (strfind(names{i}, 'noise_') > 0)])
            continue
        end
    
        try
            [a, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(case_dir, names{i}, UTCases{1, 5}), 111, 1);
        
            size(a)
        
            SLC = size(a, 3)
            for slc=1:SLC
                %if ~exist(fullfile(data_dir, [names_cases{k} '_' names{i} '_' UTCases{1, 5} '_slc' num2str(slc-1) '_PSIR.fig']))
                    h = figure('name', names{i});
                    imagescn(a(:,:,slc), [header(slc).window_center-header(slc).window_width/2 header(slc).window_center+header(slc).window_width/2], [], [12]);
               
                    saveas(h, fullfile(data_dir, [names_cases{k} '_' names{i} '_' UTCases{1, 5} '_slc' num2str(slc-1) '_PSIR.fig']), 'fig');
                    close(h)
                %end
            end            
        catch
        end
    end
end

% pick the cases

picked_dir = '/fastdata/data/db_lge/picked/'
picked_dir = '/fastdata/data/db_lge/picked_saved'

for k=1:num_cases

    k
    case_dir = fullfile(data_dir, names_cases{k})

    findAndMoveMeasDat(case_dir)
    
    UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;
    
    [names, num] = FindSubDirs(case_dir);

    UTCases{1, 4} = 'GTPrep_DB_LGE.xml';
    UTCases{1, 5} = 'res_GTPrep_DB_LGE'; 

    for i=1:num
        if any([(strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0) (strfind(names{i}, 'noise_') > 0)])
            continue
        end
    
        try
            [a, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(case_dir, names{i}, UTCases{1, 5}), 111, 1);
        
            size(a)
        
            SLC = size(a, 3)
            for slc=1:SLC
                h = figure('name', names{i});
                imagescn(a(:,:,slc), [header(slc).window_center-header(slc).window_width/2 header(slc).window_center+header(slc).window_width/2], [], [12]);
                
                reply = input('Do you want to keep this case? y/n [n]:','s');
                if isempty(reply)
                  reply = 'n';
                end
                
                if strcmp(reply, 'y')
                    saveas(h, fullfile(picked_dir, [names_cases{k} '_' names{i} '_' UTCases{1, 5} '_slc' num2str(slc-1) '_PSIR.fig']), 'fig');
                end
                close(h)
            end            
        catch
        end
    end
end

% find all picked cases

cases_names = findFILE(picked_dir, '*.fig');

N = numel(cases_names);

case_dirs = cell(N, 1);
h5names = cell(N, 1);
slcs = cell(N, 1);

for c=1:numel(cases_names)

    fig_name = cases_names{c}

    [p, fname, ext] = fileparts(fig_name);

    ind = strfind(fname, '--');

    case_dir = fname(1:ind(1)-1);
    h5name = fname(ind(1)+2:ind(2)-1);

    ind3 = strfind(fname, '_slc');
    ind4 = strfind(fname, '_PSIR');

    slc = str2num(fname(ind3+4:ind4-1));

    case_dirs{c, 1} = case_dir;
    h5names{c, 1} = h5name;
    slcs{c, 1} = slc;
end

unique_case_dirs = {}
unique_h5names = {}
for c=1:numel(cases_names)
    already_exists = 0;
    for k=1:numel(unique_h5names)
        if strcmp(h5names{c, 1}, unique_h5names{k})
            already_exists = 1;
            break;
        end
    end

    if ~already_exists
        unique_h5names = [unique_h5names {h5names{c, 1}}];
        unique_case_dirs = [unique_case_dirs {case_dirs{c, 1}}];
    end
end

numel(unique_h5names)

for c=1:numel(unique_h5names)

    disp([num2str(c) ' - ' fullfile(data_dir, unique_case_dirs{c}, unique_h5names{c})])

    UTCases = set_up_UT_cases_FreeMax_AI_Denoising_v2;
        
    UTCases{1, 1} = fullfile(data_dir, unique_case_dirs{c});
    UTCases{1, 2} = unique_h5names{c};
    
    for p=1:6
        if p==1
            UTCases{1, 4} = 'GTPrep_DB_LGE_STCNNT.xml';
            UTCases{1, 5} = 'res_GTPrep_DB_LGE_STCNNT';
        end
        if p==2
            UTCases{1, 4} = 'GTPrep_DB_LGE.xml';
            UTCases{1, 5} = 'res_GTPrep_DB_LGE'; 
        end
        if p==3
            UTCases{1, 4} = 'GTPrep_DB_LGE_half.xml';
            UTCases{1, 5} = 'res_GTPrep_DB_LGE_half'; 
        end
        if p==4
            UTCases{1, 4} = 'GTPrep_DB_LGE_STCNNT_half.xml';
            UTCases{1, 5} = 'res_GTPrep_DB_LGE_STCNNT_half'; 
        end
        if p==5
            UTCases{1, 4} = 'GTPrep_DB_LGE_STCNNT_quarter.xml';
            UTCases{1, 5} = 'res_GTPrep_DB_LGE_STCNNT_quarter'; 
        end
        if p==6
            UTCases{1, 4} = 'GTPrep_DB_LGE_quarter.xml';
            UTCases{1, 5} = 'res_GTPrep_DB_LGE_quarter'; 
        end
        
        % fnames = findFILE(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, UTCases{1, 5}), '*.hdr');
        % if (numel(fnames)>0)
        %     continue
        % end       
        
        performUTValidation(UTCases, 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput');

        debug_dir = fullfile(data_dir, unique_case_dirs{c}, unique_h5names{c}, UTCases{1, 5}, 'DebugOutput');
        [data, gmap, res] = prepare_for_stcnnt_inference_LGE(debug_dir); 
        % if ~isempty(res)
        %     SLC = size(data, 5)
        %     if SLC==1
        %         h = figure('Name', unique_h5names{c}); imagescn(cat(4, data, res), [], [SLC 4], [12], 3);
        %     else
        %         h = figure('Name',unique_h5names{c}); imagescn(cat(4, data, res), [], [4 SLC], [12], 3);
        %     end
        %     % saveas(h, fullfile(data_dir, [names{i} '_' UTCases{1, 5} '.fig']));
        %     close(h)
        % end
    end
end

%% make the mat files

mat_dir = '/fastdata/data/db_lge/picked_saved_2/'
mkdir(mat_dir)

moco_key_frame = 3 + 1;

for c=46:numel(cases_names)

    disp([num2str(c) ' - ' fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1})])

    [psir_ai_no_ave, header_psir_ai_no_ave, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_STCNNT'), 107, 1);
    [psir_ai_half_no_ave, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_STCNNT_half'), 107, 1);
    [psir_ai_quarter_no_ave, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_STCNNT_quarter'), 107, 1);

    [psir_no_ave, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE'), 107, 1);
    [psir_half_no_ave, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_half'), 107, 1);
    [psir_quarter_no_ave, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_quarter'), 107, 1);

    size(psir_ai_no_ave)

    [psir_ai, header_psir_ai, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_STCNNT'), 111, 1);
    [psir_ai_half, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_STCNNT_half'), 111, 1);
    [psir_ai_quarter, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_STCNNT_quarter'), 111, 1);

    [psir, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE'), 111, 1);
    [psir_half, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_half'), 111, 1);
    [psir_quarter, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_quarter'), 111, 1);

    size(psir_ai)

    slc = slcs(c);
    slc = slc{1} + 1;

    %im = cat(3, psir(:,:,slc), psir_ai(:,:,slc), psir_half(:,:,slc), psir_ai_half(:,:,slc), psir_quarter(:,:,slc), psir_ai_quarter(:,:,slc));
    
    im = cat(3, psir(:,:,slc), psir_half(:,:,slc), psir_quarter(:,:,slc), ...
        psir_ai(:,:,slc), psir_ai_half(:,:,slc), psir_ai_quarter(:,:,slc), ...
        psir_ai_no_ave(:,:,moco_key_frame, slc), psir_ai_half_no_ave(:,:,moco_key_frame, slc), psir_ai_quarter_no_ave(:,:,min([moco_key_frame, size(psir_ai_quarter_no_ave, 3)]), slc) ...
        );

    min_psir_ai = header_psir_ai(slc).window_center-header_psir_ai(slc).window_width/2;
    max_psir_ai = header_psir_ai(slc).window_center+header_psir_ai(slc).window_width/2;
    min_psir_ai_no_ave = header_psir_ai_no_ave(slc).window_center-header_psir_ai_no_ave(slc).window_width/2;
    max_psir_ai_no_ave = header_psir_ai_no_ave(slc).window_center+header_psir_ai_no_ave(slc).window_width/2;

    scales = [min_psir_ai max_psir_ai; min_psir_ai max_psir_ai; min_psir_ai max_psir_ai;  min_psir_ai max_psir_ai; min_psir_ai max_psir_ai; min_psir_ai max_psir_ai; ...
        min_psir_ai_no_ave max_psir_ai_no_ave;  min_psir_ai_no_ave max_psir_ai_no_ave;  min_psir_ai_no_ave max_psir_ai_no_ave];


    h = figure('name', [h5names{c, 1} '-slc' num2str(slc-1)]);
    imagescn(permute(im, [2 1 3]), [scales], [3 3], [8]);
    
    debug_dir = fullfile(data_dir, case_dirs{c, 1}, h5names{c, 1}, 'res_GTPrep_DB_LGE_STCNNT', 'DebugOutput');
    [data, gmap, res] = prepare_for_stcnnt_inference_LGE(debug_dir); 

    psir = squeeze(data(:,:,:,1,:) .* conj(data(:,:,:,2,:)./abs(data(:,:,:,2,:))));
    psir_ai = squeeze(res(:,:,:,1,:) .* conj(res(:,:,:,2,:)./abs(res(:,:,:,2,:))));

    h2 = figure('Name',  h5names{c, 1}); imagescn(real(cat(4, psir(:,:,:,slc), psir_ai(:,:,:,slc))), [-10 10], [1 2], [16], 3);

    % reply = input('Do you want to keep this case? y/n [n]:','s');
    % if isempty(reply)
    %   reply = 'n';
    % end
    
    reply = 'y';
    if strcmp(reply, 'y')
        saveas(h, fullfile(mat_dir, [case_dirs{c, 1} '--' h5names{c, 1} '--' UTCases{1, 5} '_slc' num2str(slc-1) '_PSIR.fig']), 'fig');
        saveas(h2, fullfile(mat_dir, [case_dirs{c, 1} '--' h5names{c, 1} '--' UTCases{1, 5} '_slc' num2str(slc-1) '_PSIR_no_moco_ave.fig']), 'fig');
        save(fullfile(mat_dir, [case_dirs{c, 1} '--' h5names{c, 1} '--' UTCases{1, 5} '_slc' num2str(slc-1) '_PSIR.mat']), 'im', 'psir', 'psir_ai');
    end

    close(h)
    close(h2)
end

%% rtcine

data_dir = '/data/raw_data/rtcine_R5/'
[cases, n_cases] = FindSubDirs(data_dir)

res_dir = '/data/raw_data/rtcine_R5_res'
mkdir(res_dir)

UTCases = set_up_UT_cases_RTCine;

for n=39:n_cases
    case_dir = fullfile(data_dir, cases{n});

    findAndMoveMeasDat(case_dir)
    
    [names, num] = FindSubDirs(case_dir)

    res_case_dir = fullfile(res_dir, cases{n});
    mkdir(res_case_dir)

    for i=1:num
        UTCases{1, 1} = case_dir;
        UTCases{1, 2} = names{i};
       
        UTCases{1, 4} = 'Generic_RTCine_STCNNT.xml';
        UTCases{1, 5} = 'res_Generic_RTCine_STCNNT';
    
         if any([(strfind(names{i}, 'ISMRMRD_Noise_dependency') > 0) (strfind(names{i}, 'noise_') > 0)])
                UTCases{1, 4} = 'default_measurement_dependencies.xml';
         end
    
        performUTValidation(UTCases(1,:), 0, 0, 'localhost', '9002', 1, 1, 0, 0, 0, [], [], [], '/home/xueh/Debug/DebugOutput')
        names{i}
    
        try

            [data, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(case_dir, names{i}, UTCases{1, 5}), 1, 1);
            [res, header, acq_time, physio_time, endo_pt, epi_pt, user_int] = readGTPlusExportImageSeries_Squeeze(fullfile(case_dir, names{i}, UTCases{1, 5}), 2101, 1);

            dset = ismrmrd.Dataset(fullfile(case_dir, names{i}, [names{i} '.h5']), 'dataset');
            hdr = ismrmrd.xml.deserialize(dset.readxml);
            dset.close()
    
            if ndims(data)==4
                data = permute(data, [1 2 4 3]);
                res = permute(res, [1 2 4 3]);
            end

            h = figure('Name', [names{i} '_' hdr.measurementInformation.protocolName '_R' num2str(hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1) '_RO' num2str(hdr.encoding.reconSpace.matrixSize.x)]); 
            if size(data, 4) <= 2
                imagescn(cat(4, data, res), [0 3*abs(mean(data(:)))], [size(data, 4), 2], [12], 3);
            else
                imagescn(cat(4, data, res), [0 3*abs(mean(data(:)))], [2 size(data, 4)], [12], 3);
            end
            saveas(h, fullfile(res_case_dir, [names{i} '_' regexprep(hdr.measurementInformation.protocolName, '\s+', '_') '_R' num2str(hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1) '_RO' num2str(hdr.encoding.reconSpace.matrixSize.x) '.fig']));

            save(fullfile(res_case_dir, [names{i} '.mat']), 'data', 'res', 'hdr');

            [dicom_cases, n_dicom] = FindSubDirs('/tmp/gadgetron_data');
            for tt=1:n_dicom
                movefile(fullfile('/tmp/gadgetron_data', dicom_cases{tt}), res_case_dir);
            end

            close(h);
            
            % debug_dir = fullfile(data_dir, names{i}, UTCases{1, 5}, 'DebugOutput');
            % [data, gmap, res] = prepare_for_stcnnt_inference_cine(debug_dir); 
            % h = figure('Name', [names{i} '_' hdr.measurementInformation.protocolName '_R' num2str(hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1) '_RO' num2str(hdr.encoding.reconSpace.matrixSize.x)]); 
            % if size(data, 4) <= 2
            %     imagescn(cat(4, data, res), [0 3*abs(mean(data(:)))], [size(data, 4), 2], [12], 3);
            % else
            %     imagescn(cat(4, data, res), [0 3*abs(mean(data(:)))], [2 size(data, 4)], [12], 3);
            % end
            % saveas(h, fullfile(data_dir, [names{i} '_' hdr.measurementInformation.protocolName '_R' num2str(hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1) '_RO' num2str(hdr.encoding.reconSpace.matrixSize.x) '.fig']));
        catch
        end
    end
end

%% fast MRI

dat_dir = '/data/raw_data/fastmri/knee/multicoil_train/'
names = findFILE(dat_dir, '*.h5');

res_dir = '/data/raw_data/fastmri/multicoil_train_for_model/'

R = [2 3 4]

refE1 = 32;

for k=1:numel(names)
   
    k
    [d, fname, ext] = fileparts(names{k})
    
    h = h5read(names{k}, '/ismrmrd_header');
    kspace = h5read(names{k}, '/kspace');
    
    kspace = complex(kspace.r, kspace.i);
    size(kspace)
    
    kspace = permute(kspace, [2 1 3 4]);
    RO = size(kspace, 1)
    E1 = size(kspace, 2)
    CHA = size(kspace, 3)
    SLC = size(kspace, 4)
    
    ind = find(abs(kspace(E1/2, :, 1,1))>0);
    
    nn_RO = 32
    nn_E1 = 8
    
    rx_dwelltime_data = 5.0;
    rx_dwelltime_prescan_noise = 5.0;
    
    nn1 = kspace(1:nn_RO, ind(1):ind(1)+nn_E1-1, :, :);
    size(nn1)
    
    nn2 = kspace(RO-nn_RO+1:RO, ind(1):ind(1)+nn_E1-1, :, :);
    size(nn2)
    nn3 = kspace(1:nn_RO, ind(end)-nn_E1+1:ind(end), :, :);
    size(nn3)
    nn4 = kspace(RO-nn_RO+1:RO, ind(end)-nn_E1+1:ind(end), :, :);
    size(nn4)
    
    nn = cat(2, nn1, nn2, nn3, nn4);
    size(nn)
    
    nn = permute(nn, [3 1 2 4]);
    size(nn)
    
    nn = reshape(nn, [CHA nn_RO 4*nn_E1*SLC]);
    size(nn)
    
    noisePrewhiteningMatrix = calculateNoisePrewhitener2(nn, rx_dwelltime_data, rx_dwelltime_prescan_noise);
    
    inputdata = permute(kspace, [3 1 2 4]);
    outputdata = applyNoisePrewhitener(inputdata, noisePrewhiteningMatrix);
    size(outputdata)
    kspace2 = permute(outputdata, [2 3 1 4]);
    size(kspace2)
    
    complex_im = ifft2c(kspace2);
    sos_im = sum(complex_im.*conj(complex_im), 3);
    size(sos_im)
    
    %figure; imagescn(abs(sos_im), [], [], [], 4);
    
    coilMap = Matlab_gt_compute_coil_map(reshape(complex_im, [RO, E1, 1, CHA, SLC]), 'ISMRMRD_SOUHEIL', 7, 5, 15, 1e-6);
    size(coilMap)
    
    %figure; imagescn(abs(squeeze(a.coilMap)), [], [], [], 4);
    
    dst_dir = fullfile(res_dir, fname);
    mkdir(dst_dir)
   
    coilMap = squeeze(coilMap);
    im = sum(complex_im .* conj(coilMap), 3);
    im = squeeze(im);

    fRO = Matlab_gt_generate_kspace_filter_symmetric(RO, 'Gaussian', 1.2, 12);
    fE1 = Matlab_gt_generate_kspace_filter_symmetric(E1, 'Gaussian', 1.2, 12);
    
    filteredData = performKSpaceFilter2D(fft2c(im), fRO, fE1);
    im_fil = ifft2c(filteredData);
           
    for r=1:numel(R)
        kRO = 5;
        kE1 = 4;
        accelFactor = R(r);
        fitItself = 1;
        thres = 5e-5;
        
        data = kspace2;
        data(:) = 0;
        
        data(:,1:accelFactor:end,:,:) = kspace2(:,1:accelFactor:end,:,:);
        
        aliasedIm = ifft2c(data);
        
        snr_scaling_ratio = sqrt(numel(find(abs(kspace2)>0))/numel(find(abs(data)>0)))

        gmaps = [];
        unmixing = [];
        recon = [];
        for slc=1:SLC
            slc
            acsSrc = double(kspace2(:,E1/2-refE1/2:E1/2+refE1/2,:,slc));
            acsDst = acsSrc;
            [ker, convKer] = Matlab_gt_grappa_2d_calibrate(acsSrc, acsDst, kRO, kE1, accelFactor, fitItself, thres);
            kIm = Matlab_gt_grappa_2d_compute_image_domain_kernel(convKer, RO, E1);
            [unmixCoeff, gFactor] = Matlab_gt_grappa_2d_compute_unmxing_coeff(kIm, coilMap(:,:,:,slc), accelFactor);
        
            gmaps(:,:,slc) = gFactor;
            unmixing(:,:,:,slc) = unmixCoeff;

            a_aliased = aliasedIm(:,:,:,slc);
            de_aliased = sum(repmat(a_aliased, [1 1 1 CHA]).*kIm, 3);
            de_aliased = squeeze(fft2c(de_aliased));
            de_aliased(:,E1/2-refE1/2:E1/2+refE1/2,:) = acsSrc;

            recon(:,:,slc) = sum(ifft2c(de_aliased).*conj(coilMap(:,:,:,slc)), 3);
        end
        
        recon = recon * snr_scaling_ratio;
        % save(fullfile(dst_dir, ['res_R' num2str(accelFactor)]), 'recon', 'gmaps')

        writeNPY(single(gmaps), fullfile(dst_dir, ['gmap_R' num2str(accelFactor) '.npy']));
        writeNPY(single(real(recon)), fullfile(dst_dir, ['input_R' num2str(accelFactor) '_real.npy']));
        writeNPY(single(imag(recon)), fullfile(dst_dir, ['input_R' num2str(accelFactor) '_imag.npy']));
    end

    writeNPY(single(real(im)), fullfile(dst_dir, 'im_real.npy'));
    writeNPY(single(imag(im)), fullfile(dst_dir, 'im_imag.npy'));
end


cd /data/raw_data/fastmri/multicoil_train_for_model/file1000010/

im = complex(readNPY('im_real'), readNPY('im_imag'));
a = complex(readNPY('input_R2_real'), readNPY('input_R2_imag'));
b = complex(readNPY('input_R3_real'), readNPY('input_R3_imag'));
c = complex(readNPY('input_R4_real'), readNPY('input_R4_imag'));
size(a)

figure; imagescn(abs(cat(4, im, a, b, c)), [], [1 4], [12], 3);

g2 = readNPY('gmap_R2');
g3 = readNPY('gmap_R3');
g4 = readNPY('gmap_R4');

figure; imagescn(cat(4, g2, g3, g4), [0 6], [1 3], [12], 3);colormap('jet')

data_dir='/data/raw_data/fastmri/multicoil_train_for_model/file1000010/'

[input, res, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, 'model_res_70_R1'), 0);
[input_R2, res_R2, gmap] = load_results_stcnnt_inference_perf(fullfile(data_dir, 'model_res_70_R2'), 0);

figure; imagescn(abs(cat(4, input, res, res_R2)), [], [1 3], [12], 3);

%% run the cases

[cases, ncases] = FindSubDirs(res_dir)

run_str = '_R1'
scaling_factor = 7.5

for g=1:4

    gpu_str = num2str(g-1)
    check_processed = 0
    a_model = '/data/raw_data/retro_cine_3T/models/ws8_8_16_ps4_4_2_patch_size4_B1_2xG8_epoch_89_neuro.pth'
    a_model = '/data/models/OR_T1L1G1T1L1G1_64_asp_1_dp_0.1_epoch_denoising_no_dealiasing_layer_nn_128.0_16.0_8xV100/kings01_hrnet__8xG16/checkpoints/recon-STCNNT_HRNET_kings01_hrnet__8xG16_8x32G16-V100_NN_128.0_16.0_C-64_complex_residual/checkpoint_epoch_159.pth'
    batch_size = '1'
    
    code_dir = '/home/xueh/mrprogs/imagingfm_BTCHW'
    a_run_name = fullfile(res_dir, ['run' run_str '_' gpu_str '.sh']);
    
    f = fopen(a_run_name, "w")
    fprintf(f, '%s\n', "#!/bin/bash");
    fprintf(f, '%s\n', strjoin(["export CUDA_VISIBLE_DEVICES=", gpu_str], ""));
    fprintf(f, '%s\n', "export DISABLE_FLOAT16_INFERENCE=True");
    fprintf(f, '%s\n', "export RES_DIR=model_res");
    fprintf(f, '%s\n', 'echo "------------------------------------------------------------------------"');
    
    for kk=g:4:ncases
        
        a_case = fullfile(res_dir, cases{kk})
        a_output_case_dir = fullfile(a_case, ['model_res' run_str])
    
        if strcmp(run_str, '_R1')
            input_name = ['im'];
        else
            input_name = ['input' run_str];
        end
        gmap_name = ['gmap' run_str];
    
        code_str = strjoin([code_dir "/projects/mri_imaging/inference/run_inference.py --check_processed " num2str(check_processed)], "");
        cmd_str = ["python3 " code_str " --input_dir "  a_case " --output_dir " a_output_case_dir " --scaling_factor " num2str(scaling_factor) " --im_scaling 1.0 --gmap_scaling 1.0 --input_fname " input_name " --gmap_fname " gmap_name " --saved_model_path " a_model " --batch_size " batch_size];
        cmd_str = strjoin(cmd_str, " ")
    
        fprintf(f, cmd_str);
        fprintf(f, '\n#############################################\n\n');
    end    
    
    fprintf(f, '\n echo "All run completed ... "\n\n');
    fclose(f);
    
    a_run_name
end

case_dir = '/data/raw_data/fastmri/multicoil_train_for_model/file100001'

im = complex(readNPY(fullfile(case_dir, 'im_real')), readNPY(fullfile(case_dir, 'im_imag.npy')));

[input1, res1, gmap] = load_results_stcnnt_inference_perf(fullfile(case_dir, 'model_res_R3'), 0);
[input2, res2, gmap] = load_results_stcnnt_inference_perf(fullfile(case_dir, 'model_res_R4'), 0);
[input3, res3, gmap] = load_results_stcnnt_inference_perf(fullfile(case_dir, 'model_res_R1'), 0);

figure; imagescn(cat(4, abs(im), res3), [], [1 2], [12], 3);

figure; imagescn(cat(4, abs(im), input1, res1, input2, res2), [], [1 5], [12], 3);

case_dir = '/data/raw_data/fastmri/multicoil_train_for_model/file1000013/'
case_dir = '/data/raw_data/fastmri/multicoil_train_for_model/file1000010/'
case_dir = '/data/raw_data/fastmri/multicoil_train_for_model/file1002262/'

im = complex(readNPY(fullfile(case_dir, 'im_real')), readNPY(fullfile(case_dir, 'im_imag.npy')));
[input3, res3, gmap] = load_results_stcnnt_inference_perf(fullfile(case_dir, 'model_res_R1'), 0);
figure; imagescn(cat(4, abs(im), res3), [], [1 2], [12], 3);

%% create h5 for training
% python3 ./projects/mri_imaging/data/create_hdf5_fastmri_dataset.py --output /data/raw_data/fastmri/multicoil_train_for_model/knee.h5 /data/raw_data/fastmri/multicoil_train_for_model --gmap_scaling 1 --im_scaling 1 --input_fname model_res_R1/output