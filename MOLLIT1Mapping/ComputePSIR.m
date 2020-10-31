function [PSIR, im_phase_removed] = ComputePSIR(im)
% compute PSIR recon
% [PSIR, im_phase_removed] = ComputePSIR(im)

backgroundPhs = im;
backgroundPhs = backgroundPhs ./ (abs(backgroundPhs)+eps);

im_phase_removed = im .* conj(backgroundPhs);
PSIR = real(im_phase_removed);