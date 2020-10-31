
function Kmeans_InitialCentres = ReadKmeans_InitialCentres(Kmeansfile, classnumber)
% read kmeans file

if ( classnumber == 4 )
    size = 4;
else
    size = 5;
end

fid = fopen(Kmeansfile, 'r');
[Kmeans_InitialCentres,count] = fscanf(fid,'%f',size);
fclose(fid);
return