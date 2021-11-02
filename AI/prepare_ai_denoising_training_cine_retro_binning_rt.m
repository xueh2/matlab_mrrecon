function prepare_ai_denoising_training_cine_retro_binning_rt(dataDir, resDir, resDir_local, aiDir, files_all, im_series_num, gfactor_series_num, check_processed)
% prepare_ai_denoising_training_cine_retro_binning_rt(dataDir, resDir, resDir_local, aiDir, files_all, im_series_num, gfactor_series_num, check_processed)

mkdir(aiDir)

pic_dir = fullfile(aiDir, 'jpg_pics');
mkdir(pic_dir)

pic_data_dir = fullfile(pic_dir, 'data');
mkdir(pic_data_dir)
pic_gmap_dir = fullfile(pic_dir, 'gmap');
mkdir(pic_gmap_dir)

accelFactor = [2 3 4 5];

for n = 1:size(files_all,1)
    
    fname = files_all{n, 1};
    if(iscell(fname))
        fname = fname{1};
    end
    [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(fname);
    
    h5_name = fullfile(dataDir, study_dates, [fname '.h5']);
    
    case_dir = fullfile(resDir, study_dates, fname);   
    case_dir_debug = case_dir;
    if(~isempty(resDir_local))
        case_dir_debug = fullfile(resDir_local, study_dates, fname);
        if(~exist(case_dir_debug))
            case_dir_debug = case_dir;
        end
    end
    dst_dir = fullfile(aiDir, study_dates, fname);
    
%     [configName, scannerID, patientID, studyID, measurementID, study_dates, study_year, study_month, study_day, study_time] = parseSavedISMRMRD(files_all{n});
%     
%     case_dir = fullfile(resDir, study_dates, files_all{n});    
%     dst_dir = fullfile(aiDir, study_dates, files_all{n});
    
    if(exist(case_dir))        
    
        disp(['--> process ' num2str(n) ' out of ' num2str(size(files_all,1)) ' - ' files_all{n}]);
        
        if(check_processed & exist(fullfile(pic_gmap_dir, [fname '.jpg'])))
            continue;
        end
        
        debug_dir = fullfile(case_dir_debug, 'DebugOutput');
        dset = ismrmrd.Dataset(h5_name, 'dataset');
        hdr = ismrmrd.xml.deserialize(dset.readxml);
        dset.close();

        gmap_slices = [];

        SLC = hdr.encoding(1).encodingLimits.slice.maximum + 1;
                
%         if(exist(fullfile(dst_dir, 'im_real.npy')))
%             im_real = readNPY(fullfile(dst_dir, 'im_real.npy'));
%             im_imag = readNPY(fullfile(dst_dir, 'im_imag.npy'));
%             im = complex(im_real, im_imag);
%             gfactor = readNPY(fullfile(dst_dir, 'gfactor.npy'));
%             im = abs(im);
%         else
%             try
%                 [im, im_header, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, im_series_num, 1);
%                 [gfactor, gmap_header, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, gfactor_series_num, 1);
%             catch
%                 continue;                
%             end
%             im = squeeze(im);
%             gfactor = squeeze(gfactor);
%             gfactor = abs(gfactor);
% 
%             mkdir(dst_dir);
% 
%             if(isreal(im))
%                 writeNPY(single(im), fullfile(dst_dir, 'im.npy'));
%             else
%                 writeNPY(single(real(im)), fullfile(dst_dir, 'im_real.npy'));
%                 writeNPY(single(imag(im)), fullfile(dst_dir, 'im_imag.npy'));
%                 writeNPY(single(abs(complex(im_real, im_imag))), fullfile(dst_dir, 'im.npy'));
%            end
%             writeNPY(single(gfactor), fullfile(dst_dir, 'gfactor.npy'));
%             save(fullfile(dst_dir, 'headers.mat'), 'im_header', 'gmap_header');
%             
%             im = abs(im);
%         end
%         
%         im = uint8(255 * im / (0.5*max(im(:))) );
%         
%         if(numel(size(gfactor))==2)
%             gfactor = repmat(gfactor, [1 1 size(im, 3)]);
%         end
%         
%         if(numel(size(gfactor))==3 & numel(size(im))==4)
%             gfactor = repmat(gfactor, [1 1 1 size(im, 4)]);
%         end
%         
%         if(numel(size(im))==3)
%             h = figure; imagescn(cat(4, im, gfactor), [], [1 2], [12], 3);
%         else
%             SLC = size(im, 3);
%             h = figure; imagescn(cat(3, im, gfactor), [], [2 SLC], [24], 4);
%         end
%         saveas(h, fullfile(pic_dir, [files_all{n} '.jpg']), 'jpg');
%         closeall

        try
            [im, im_header, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, im_series_num, 1);
            [gfactor, gmap_header, acq_time, physio_time] = readGTPlusExportImageSeries_Squeeze(case_dir, gfactor_series_num, 1);
        catch
            continue;                
        end
        im = squeeze(im);
        gfactor = squeeze(gfactor);

        mkdir(dst_dir);

        if(isreal(im))
            writeNPY(single(im), fullfile(dst_dir, 'im.npy'));
        else
            writeNPY(single(real(im)), fullfile(dst_dir, 'im_real.npy'));
            writeNPY(single(imag(im)), fullfile(dst_dir, 'im_imag.npy'));
            writeNPY(single(abs(im)), fullfile(dst_dir, 'im.npy'));
        end
        writeNPY(single(gfactor), fullfile(dst_dir, 'gfactor.npy'));
        save(fullfile(dst_dir, 'headers.mat'), 'im_header', 'gmap_header');
        
        try
            gmap_slices = [];
            acs_src = readGTPlusExportData(fullfile(debug_dir, ['ref_calib_encoding_0_' num2str(1)]));
            acs_dst = readGTPlusExportData(fullfile(debug_dir, ['ref_calib_dst_encoding_0_' num2str(1)]));

            acs_dst_coil_map = readGTPlusExportData(fullfile(debug_dir, ['ref_coil_map_dst_encoding_0_' num2str(1)]));

            acs_src = squeeze(acs_src);
            acs_dst = squeeze(acs_dst);
            acs_dst_coil_map = squeeze(acs_dst_coil_map);

            RO = size(acs_dst_coil_map, 1);
            E1 = size(acs_dst_coil_map, 2);
            srcCHA = size(acs_src, 3);
            dstCHA = size(acs_dst, 3);

            c_im = ifft2c(acs_dst_coil_map);
                
            for slc=1:SLC
                disp(['--> process slice ' num2str(slc)]);
                coil_map = Matlab_gt_compute_coil_map(single(reshape(c_im(:,:,:,slc), [RO E1 1 dstCHA 1])), 'ISMRMRD_SOUHEIL', 7, 5, 5, 1e-3);
                coil_map = squeeze(coil_map);

                clear gmap
                for af=1:numel(accelFactor)
                    kRO = 5;
                    kE1 = 4;
                    fitItself = 1;
                    thres = 5e-4;
                    [ker, convKer] = Matlab_gt_grappa_2d_calibrate(double(acs_src(:,:,:,slc)), double(acs_dst(:,:,:,slc)), kRO, kE1, accelFactor(af), fitItself, thres);

                    kIm = Matlab_gt_grappa_2d_compute_image_domain_kernel(double(convKer), RO, E1);

                %     figure; imagescn(abs(kIm), [], [], [], 3)

                    [unmixing, gmap(:,:,af)] = Matlab_gt_grappa_2d_compute_unmxing_coeff(kIm, coil_map, accelFactor(af));
                end

                gmap = Matlab_gt_resize_2D_image(gmap, size(im,1), size(im,2), 5);

                writeNPY(single(gmap), fullfile(dst_dir, ['gmap_slc_' num2str(slc) '.npy']));

                gmap_slices(:,:,:,slc) = gmap;
            end
        catch
            gmap_slices = [];
            for slc=1:SLC
                disp(['--> process slice ' num2str(slc)]);
                
                acs_src = readGTPlusExportData(fullfile(debug_dir, ['acsSrc_n_0s_' num2str(slc-1)]));
                acs_dst = readGTPlusExportData(fullfile(debug_dir, ['acsDst_n_0s_' num2str(slc-1)]));

%                 acs_dst_coil_map = readGTPlusExportData(fullfile(debug_dir, ['coilMap_n_0s_' num2str(slc-1)]));

                acs_src = squeeze(acs_src);
                acs_dst = squeeze(acs_dst);
%                 acs_dst_coil_map = squeeze(acs_dst_coil_map);
% 
%                 RO = size(acs_dst_coil_map, 1);
%                 E1 = size(acs_dst_coil_map, 2);
%                 srcCHA = size(acs_src, 3);
%                 dstCHA = size(acs_dst, 3);
% 
%                 c_im = ifft2c(acs_dst_coil_map);
%             
%                 coil_map = Matlab_gt_compute_coil_map(single(reshape(c_im, [RO E1 1 dstCHA 1])), 'ISMRMRD_SOUHEIL', 7, 5, 5, 1e-3);
%                 coil_map = squeeze(coil_map);
% 
                coil_map = readGTPlusExportData(fullfile(debug_dir, ['coilMap_n_0s_' num2str(slc-1)]));
                coil_map = squeeze(coil_map);
                
                RO = size(coil_map, 1);
                E1 = size(coil_map, 2);
                srcCHA = size(acs_src, 3);
                dstCHA = size(acs_dst, 3);
                
                clear gmap
                for af=1:numel(accelFactor)
                    kRO = 5;
                    kE1 = 4;
                    fitItself = 1;
                    thres = 5e-4;
                    
                    [ker, convKer] = Matlab_gt_grappa_2d_calibrate(double(acs_src), double(acs_dst), kRO, kE1, accelFactor(af), fitItself, thres);

                    kIm = Matlab_gt_grappa_2d_compute_image_domain_kernel(double(convKer), RO, E1);

                %     figure; imagescn(abs(kIm), [], [], [], 3)

                    [unmixing, gmap(:,:,af)] = Matlab_gt_grappa_2d_compute_unmxing_coeff(kIm, coil_map, accelFactor(af));
                end

                gmap = Matlab_gt_resize_2D_image(gmap, size(im,1), size(im,2), 5);

                writeNPY(single(gmap), fullfile(dst_dir, ['gmap_slc_' num2str(slc) '.npy']));

                gmap_slices(:,:,:,slc) = gmap;
            end
        end

        im = abs(im);
        gfactor = abs(gfactor);

        im = uint8(255 * im / (0.5*max(im(:))) );

        if(numel(size(gfactor))==2)
            gfactor = repmat(gfactor, [1 1 size(im, 3)]);
        end

        if(numel(size(gfactor))==3 & numel(size(im))==4)
            gfactor = repmat(gfactor, [1 1 1 size(im, 4)]);
        end

        if(numel(size(im))==3)
            h = figure; imagescn(cat(4, im, gfactor), [], [1 2], [12], 3);
        else
            SLC = size(im, 3);
            h = figure; imagescn(cat(3, im, gfactor), [], [2 SLC], [24], 4);
        end
        saveas(h, fullfile(pic_data_dir, [fname '.jpg']), 'jpg');

        try
            RO = size(im, 1);
            E1 = size(im, 2);
            im_slices = im(:,:,:,1);
            d = zeros(RO, E1, numel(accelFactor)+1,SLC);
            if(SLC>1)
                d(:,:,1,:) = im_slices(:,:,:);
            else
                d(:,:,1,:) = im_slices(:,:,1);
            end
            d(:,:,2:end,:) = gmap_slices;
            if(SLC<3)
                h = figure; imagescn(d, [], [SLC, 5], [14]);
            else
                h = figure; imagescn(d, [], [ceil(SLC/2), 10], [14]);
            end
            saveas(h, fullfile(pic_gmap_dir, [fname '.jpg']), 'jpg');
        catch
        end
        closeall
    end    
end