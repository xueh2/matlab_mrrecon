function h = ComputeDicomCoordFromGtOffline(position, read_dir, phase_dir, field_of_view, matrix_size, RO, E1)
% h = ComputeDicomCoordFromGtOffline(position, read_dir, phase_dir, field_of_view, matrix_size, RO, E1)

% compute dicom position
posVec = ComputeDicomPosVecFromICE(position, read_dir, phase_dir, field_of_view(1)/matrix_size(1), field_of_view(2)/matrix_size(2), E1, RO);

h.PixelSpacing(1,1) = field_of_view(1)/matrix_size(1);
h.PixelSpacing(2,1) = field_of_view(2)/matrix_size(2);

h.SliceThickness = field_of_view(3)/matrix_size(3);

h.ImagePositionPatient = posVec(:);

h.ImageOrientationPatient(1:3,1) = read_dir';
h.ImageOrientationPatient(4:6,1) = phase_dir';

h.Columns = RO;
h.Rows = E1;