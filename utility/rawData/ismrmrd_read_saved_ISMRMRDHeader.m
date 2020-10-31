
function h = ismrmrd_read_saved_ISMRMRDHeader(header_file)
% h = ismrmrd_read_saved_ISMRMRDHeader(header_file)

fid = fopen(header_file, 'r');

h.version = fread(fid, 1, 'uint16');
h.data_type = fread(fid, 1, 'uint16');
h.flags = fread(fid, 1, 'uint64');
h.measurement_uid = fread(fid, 1, 'uint32');
h.matrix_size = fread(fid, 3, 'uint16');
h.field_of_view = fread(fid, 3, 'single');
h.channels = fread(fid, 1, 'uint16');
h.position = fread(fid, 3, 'single');
h.read_dir = fread(fid, 3, 'single');
h.phase_dir = fread(fid, 3, 'single');
h.slice_dir = fread(fid, 3, 'single');
h.patient_table_position = fread(fid, 3, 'single');
h.average = fread(fid, 1, 'uint16');
h.slice = fread(fid, 1, 'uint16');
h.contrast = fread(fid, 1, 'uint16');
h.phase = fread(fid, 1, 'uint16');
h.repetition = fread(fid, 1, 'uint16');
h.set = fread(fid, 1, 'uint16');
h.acquisition_time_stamp = fread(fid, 1, 'uint32');
h.physiology_time_stamp = fread(fid, 3, 'uint32');
h.image_type = fread(fid, 1, 'uint16');
h.image_index = fread(fid, 1, 'uint16');
h.image_series_index = fread(fid, 1, 'uint16');
h.user_int = fread(fid, 8, 'uint32');
h.user_float = fread(fid, 8, 'single');
h.attribute_string_len = fread(fid, 8, 'uint32');

h.PatientPosition = h.position;
h.FOV = h.field_of_view;
h.window_center = -1;
h.window_width = -1;
h.TI = 0;
h.slice_location = dot(h.slice_dir, h.position);
     
% v = h.phase_dir;
% h.phase_dir = h.read_dir;
% h.read_dir = v
% h.slice_dir = -h.slice_dir;

fclose(fid);


