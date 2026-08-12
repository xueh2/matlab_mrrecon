function [header, data, ref_data] = load_mrd_kspace(fname)
% [header, data, ref_data] = load_mrd_kspace(fname)

r = mrd.binary.MrdReader(fname)
header = r.read_header();

matrix_size = header.encoding.encoded_space.matrix_size;

data = [];
ref_data = [];
i = 1;
while r.has_data()
    item = r.read_data();
    acq = item.value;

    if acq.head.flags.has_flags(mrd.AcquisitionFlags.IS_NOISE_MEASUREMENT)
        continue
    end

    e1 = acq.head.idx.kspace_encode_step_1;
    e2 = acq.head.idx.kspace_encode_step_2;

    CHA = size(acq.data, 1);
    RO = acq.head.center_sample * 2;

    if isempty(data)
        data = zeros(CHA, matrix_size.z, matrix_size.y, RO);
    end

    if acq.head.flags.has_flags(mrd.AcquisitionFlags.IS_PARALLEL_CALIBRATION) | acq.head.flags.has_flags(mrd.AcquisitionFlags.IS_PARALLEL_CALIBRATION_AND_IMAGING)
        if isempty(ref_data)
            ref_data = zeros(CHA, matrix_size.z, matrix_size.y, acq.head.center_sample * 2);
            ref_RO = acq.head.center_sample * 2;
        end

        ref_data(:, e2+1, e1+1, ref_RO-size(acq.data,2)+1:end) = acq.data;
    end

    if ~acq.head.flags.has_flags(mrd.AcquisitionFlags.IS_PARALLEL_CALIBRATION)
        data(:, e2+1, e1+1, RO-size(acq.data,2)+1:end) = acq.data;
    end

    
    i = i + 1;
end
r.close();

disp(['load ' fname ' - ' num2str(size(data))]);