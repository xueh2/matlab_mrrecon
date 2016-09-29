function s = ismrmrd_check_readout_flags(readouts)
% check the status of a readout
% enum ISMRMRD_AcquisitionFlags {
%     ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP1               =  1,
%     ISMRMRD_ACQ_LAST_IN_ENCODE_STEP1                =  2,
%     ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP2               =  3,
%     ISMRMRD_ACQ_LAST_IN_ENCODE_STEP2                =  4,
%     ISMRMRD_ACQ_FIRST_IN_AVERAGE                    =  5,
%     ISMRMRD_ACQ_LAST_IN_AVERAGE                     =  6,
%     ISMRMRD_ACQ_FIRST_IN_SLICE                      =  7,
%     ISMRMRD_ACQ_LAST_IN_SLICE                       =  8,
%     ISMRMRD_ACQ_FIRST_IN_CONTRAST                   =  9,
%     ISMRMRD_ACQ_LAST_IN_CONTRAST                    = 10,
%     ISMRMRD_ACQ_FIRST_IN_PHASE                      = 11,
%     ISMRMRD_ACQ_LAST_IN_PHASE                       = 12,
%     ISMRMRD_ACQ_FIRST_IN_REPETITION                 = 13,
%     ISMRMRD_ACQ_LAST_IN_REPETITION                  = 14,
%     ISMRMRD_ACQ_FIRST_IN_SET                        = 15,
%     ISMRMRD_ACQ_LAST_IN_SET                         = 16,
%     ISMRMRD_ACQ_FIRST_IN_SEGMENT                    = 17,
%     ISMRMRD_ACQ_LAST_IN_SEGMENT                     = 18,
%     ISMRMRD_ACQ_IS_NOISE_MEASUREMENT                = 19,
%     ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION             = 20,
%     ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING = 21,
%     ISMRMRD_ACQ_IS_REVERSE                          = 22,
%     ISMRMRD_ACQ_IS_NAVIGATION_DATA                  = 23,
%     ISMRMRD_ACQ_IS_PHASECORR_DATA                   = 24,
%     ISMRMRD_ACQ_LAST_IN_MEASUREMENT                 = 25,
%     ISMRMRD_ACQ_IS_HPFEEDBACK_DATA                  = 26,
%     ISMRMRD_ACQ_IS_DUMMYSCAN_DATA                   = 27,
%     ISMRMRD_ACQ_IS_RTFEEDBACK_DATA                  = 28,
%     ISMRMRD_ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA   = 29,
% 
%     ISMRMRD_ACQ_USER1                               = 57,
%     ISMRMRD_ACQ_USER2                               = 58,
%     ISMRMRD_ACQ_USER3                               = 59,
%     ISMRMRD_ACQ_USER4                               = 60,
%     ISMRMRD_ACQ_USER5                               = 61,
%     ISMRMRD_ACQ_USER6                               = 62,
%     ISMRMRD_ACQ_USER7                               = 63,
%     ISMRMRD_ACQ_USER8                               = 64
% };
    
N = numel(readouts.data);
s = struct('is_noise', zeros(N, 1), 'is_ref', zeros(N, 1), 'is_ref_kspace', zeros(N, 1), 'is_reflect', zeros(N, 1), 'is_phase_corr', zeros(N, 1), 'is_navigator', zeros(N, 1), 'is_rt_feedback', zeros(N, 1), 'is_hp_feedback', zeros(N, 1), 'is_dummy', zeros(N, 1));

s.is_noise = readouts.head.flagIsSet(readouts.head.FLAGS.ACQ_IS_NOISE_MEASUREMENT);
s.is_ref = readouts.head.flagIsSet(readouts.head.FLAGS.ACQ_IS_NOISE_MEASUREMENT);

for n=1:N
    s.is_noise(n) = readouts.
end
N = readouts
s.is_noise = flags 