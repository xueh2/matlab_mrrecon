function [noise, dwelltime_us] = ismrmrd_get_noise(readouts)
% get the noise readouts from ismrmrd data

N = numel(readouts.data);

isNoise = readouts.head.flagIsSet(readouts.head.FLAGS.ACQ_IS_NOISE_MEASUREMENT);

noise = [];
dwelltime_us = [];
for n=1:N
    if(isNoise(n)==1)
        noise = [noise; readouts.data{n}];
        dwelltime_us = [dwelltime_us; readouts.head.sample_time_us(n)];
    end
end