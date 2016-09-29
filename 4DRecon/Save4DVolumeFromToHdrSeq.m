%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Format conversion to visualize Hui's reconstruction (separate analyze)
%
% 4DCardio project
%
% Xiaoguang Lu (xiaoguang.lu@siemens.com)
% Siemens Corporation, Corporate Research & Technology
% Mar. 2012

function Save4DVolumeFromToHdrSeq(volume, header, des_dir, suffix)

cid = suffix;

sizeT = header.sizeT;

fid_seq = fopen([des_dir '\' cid '.seq'], 'wt');
fprintf(fid_seq, '%d\n', sizeT);

for i=1:sizeT
    i
    vol = volume(:,:,:,i);
    % vol = vol * 4095 / 255;
    res = [header.spacingX header.spacingY header.spacingZ];
    fn = sprintf('%s_%04d', cid, i)
    fprintf(fid_seq, '%s\n', fn);
    saveHdrImg([des_dir '\' fn], vol, res, 'uint16');
end

fclose(fid_seq);

