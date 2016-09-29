function rxdwell=bw2rxdwell(bw, xres);
% function rxdwell=bw2rxdwell(bw, xres);
% 
% function to calculate the receiver dwell (rxdwell) in ns
% from the receiver bandwidth (bw) in Hz/pixel
% and readout resolution (xres=#pixels), given 2x readout oversampling

rxdwell=1/(2*xres*bw*1e-9); % Hz/pixel


return