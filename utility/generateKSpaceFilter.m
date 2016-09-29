function kspaceFilter = generateKSpaceFilter(filterType, filterStrength, len, sampledRange, kspaceCenterFE)

if ( len <= 1 )
    kspaceFilter = [];
    return;
end

if ( strcmp(filterType, 'None') == 0 )
    kspaceFilter = Matlab_ComputeKSpaceFilter(len, len/2, sampledRange(1)-1, sampledRange(2)-sampledRange(1)+1, kspaceCenterFE, ...
                            filterType, filterStrength, 0, 0);
    kspaceFilter = single(kspaceFilter);
else
    kspaceFilter = [];
end
