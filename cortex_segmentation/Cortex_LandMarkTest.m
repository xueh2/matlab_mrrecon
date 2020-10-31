
function [dists, dists_SDF] = Cortex_LandMarkTest(pts_file, surfaceVTKfile, levelsetFun, OriginalImageFile, OutputSDF_File)

if ( isempty(dir(surfaceVTKfile)) )
    disp(' generating surfaces ... ' );
    [data, header] = LoadAnalyze(levelsetFun, 'Real');
    mcubes2VTKFile(data, header, surfaceVTKfile, 0, [1 1 1], [0 0]);
end

if ( isempty(dir(OutputSDF_File)) )
    disp(' generating SDF function ... ' );

    accuracy = 'high';
    tMax_ReIntialize = 100;
    errorMax = 0.05;

    [SDF, header] = SignedDistranceTransform(levelsetFun, OutputSDF_File, accuracy, tMax_ReIntialize, errorMax);
else
    [SDF, header] = LoadAnalyze(OutputSDF_File, 'Real');
end

[pts_World, numCell, cells] = VTKFile2PointCell(pts_file); 

disp(' computing surface distance ... ' );
dists = LandmarkDst_Cortex(pts_World(:, 2:4), surfaceVTKfile, levelsetFun, OriginalImageFile);

disp(' computing SDF distance ... ' );
dists_SDF = LandmarkDst_Cortex_SDF(pts_World(:, 2:4), OutputSDF_File);

signs = sign(dists_SDF);
dists = dists .* signs;
return