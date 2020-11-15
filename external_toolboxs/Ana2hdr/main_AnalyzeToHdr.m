%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Format conversion to visualize Hui's reconstruction (separate analyze)
%
% 4DCardio project
%
% Xiaoguang Lu (xiaoguang.lu@siemens.com)
% Siemens Corporation, Corporate Research & Technology
% Mar. 2012

if 1
    % load E:\CMR4D\Recon4D\20120322\Im4D_BreathingGating_Linear_full.mat;
    ana_dir = 'E:\CMR4D\Recon4D\20120524\Analyze';
    des_dir = 'E:\CMR4D\Recon4D\20120524\seq\20120524';
    
    cid = '20120524';
    
    fList = dir([ana_dir '\*.hdr']);

    fid_seq = fopen([des_dir '\' cid '.seq'], 'wt');
    fprintf(fid_seq, '%d\n', length(fList));

    for i=1:length(fList)
        i
        id = fList(i).name
        idx1 = strfind(fList(i).name, 'PHS')+3;
        idx2 = strfind(fList(i).name, '.')-1;
        phs = str2num(fList(i).name(idx1:idx2));
        nii = load_nii([ana_dir '\' id]);
        vol = nii.img;
        vol = vol * 4095 / 255;
        res = nii.hdr.dime.pixdim(2:4);
        s = nii.hdr.dime.dim(2:4);
        if 0 % smaller region of interest
            vol = vol(80:220, 20:180, 50:150);
            s = size(vol);
        end
        clear nii;
        fn = sprintf('%s_%04d', cid, phs)
        fprintf(fid_seq, '%s\n', fn);
        saveHdrImg([des_dir '\' fn], vol, res, 'uint16');
    end
    
    fclose(fid_seq);

    return;
end
