function res = ifft2c_MrFtk(x)
fctr = size(x,1)*size(x,2);
res = sqrt(fctr)*fftshift(Matlab_PerformFFT2D(ifftshift(x),-1));
