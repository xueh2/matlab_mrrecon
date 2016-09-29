function bw=rxdwell2bw(rxdwell, xres);
% function bw=rxdwell2bw(rxdwell, xres);
% 
% function to calculate receiver bandwidth (bw) in Hz/pixel
% from the receiver dwell (rxdwell) in ns and readout resolution
% (xres=#pixels), given 2x readout oversampling

bw=1./(2*xres*rxdwell*1e-9); % Hz/pixel



return