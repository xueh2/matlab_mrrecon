
function WriteKmeans_InitialCentres(Kmeansfile, classnumber, Kmeans_InitialCentres)
% read kmeans file

if ( classnumber == 4 )
    Kmeans_InitialCentres = [Kmeans_InitialCentres(1) Kmeans_InitialCentres(2) ...
        (Kmeans_InitialCentres(3)+Kmeans_InitialCentres(4))/2 Kmeans_InitialCentres(5)];
end

fid = fopen(Kmeansfile, 'w');
fprintf(fid,'%d ',Kmeans_InitialCentres);
fclose(fid);
return