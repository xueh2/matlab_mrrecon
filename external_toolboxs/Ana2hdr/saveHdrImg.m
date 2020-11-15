function saveHdrImg(fn, vol, res, type)

if nargin < 4
    type = 'uint16';
end

fid_hdr = fopen([fn '.hdr'], 'wt');
fprintf(fid_hdr, '%d %d %d\n', size(vol,1), size(vol,2), size(vol,3));
fprintf(fid_hdr, '%f %f %f\n', res(1), res(2), res(3));
fprintf(fid_hdr, '30 300 -1024\n');
fprintf(fid_hdr, '0.00\n' );
fprintf(fid_hdr, '1.00 0.00 0.00 0.00\n');
fprintf(fid_hdr, '0.00 1.00 0.00 0.00\n');
fprintf(fid_hdr, '0.00 0.00 1.00 0.00\n');
switch type
    case 'uint16'
        fprintf(fid_hdr, '%d\n', 2);
        fprintf(fid_hdr, '0 %d\n', 2^12-1);
    case 'uint8'
        fprintf(fid_hdr, '%d\n', 2);
        fprintf(fid_hdr, '0 %d\n', 2^8-1);
end        
fclose(fid_hdr);

ff = fopen([fn '.img'], 'wb');
switch type
    case 'uint16'
        fwrite(ff, uint16(vol), 'uint16');
    case 'uint8'
        fwrite(ff, uint8(vol), 'uint8');
end
fclose(ff);

return;