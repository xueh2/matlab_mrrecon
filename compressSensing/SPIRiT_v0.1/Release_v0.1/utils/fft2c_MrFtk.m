function res = fft2c_MrFtk(x)
fctr = size(x,1)*size(x,2);
res = 1/sqrt(fctr)*fftshift(Matlab_PerformFFT2D(ifftshift(x),1));
