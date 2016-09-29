
%% save data as MMID4 loadable format

cd D:\gtuser\gt_windows_setup\ut\FLASH_Perfusion_HighRes\20150227_09h17m34s_29748

load PerfTestReady

N = numel(cin)
T = deltaT*[0:N-1]/1000;

data = zeros(N, 2);
data(:,1) = T;
data(:,2) = cin;

