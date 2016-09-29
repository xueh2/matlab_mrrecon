function f = hasGPU()

f = 0;
try
    h = gpuDeviceCount;
    if ( h > 0 )
        f = 1;
    end
catch
    f = 0;
end