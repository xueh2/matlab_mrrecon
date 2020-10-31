
disp('----------------------------------------------------');
disp('AllinOne_Reconstruction');

    csf_seg_filename = 'csf_seg.hdr';
    wm_seg_filename = 'wm_seg.hdr';
    
    post_wm_Real_filename = 'post_wm_Real.hdr';
    post_gm_Real_filename = 'post_gm_Real.hdr';
    post_csf_Real_filename = 'post_csf_Real.hdr';
    Brainfilename = 'sum.hdr';


    Cortex_Reconstruction_WholePipeline_AllRun
disp('----------------------------------------------------');


    %------------------------------------------------------------------
    filename = 'N3Brain_roi.hdr';
    [imagedata, header] = LoadAnalyze(filename, 'Grey');

    filename = 'Internal_levelset_Result.hdr';
    [Internal, header] = LoadAnalyze(filename, 'Real');

    filename = 'External_levelset_Result_TH.hdr';
    [External, header] = LoadAnalyze(filename, 'Real');
    
    xsize = header.xsize;
    ysize = header.ysize;
    zsize = header.zsize;

    xvoxelsize = header.xvoxelsize;
    yvoxelsize = header.yvoxelsize;
    zvoxelsize = header.zvoxelsize;

    % perform the level set propagation
    %--------------------------------------------------------------------------
    accuracy = 'medium';
    tMax_ReIntialize = 100;
    errorMax = 0.05;
    %--------------------------------------------------------------------------
    CortexThickness_SignedDistanceFunction
    
    

vtkfile = 'Internal_cortex.vtk';
threshold = 0;
shrink = [1 1 1];
smooth = [0.1 20];

[Internal, header] = LoadAnalyze('Internal_levelset_Result.hdr', 'Real');
mcubes2VTKFile(Internal, header, vtkfile, threshold, shrink, smooth);

vtkfile = 'External_cortex.vtk';
threshold = 0;
shrink = [1 1 1];
smooth = [0.1 20];

[Internal, header] = LoadAnalyze('External_levelset_Result_TH.hdr', 'Real');
mcubes2VTKFile(Internal, header, vtkfile, threshold, shrink, smooth);

vtkfile = 'Middle_cortex.vtk';
threshold = 0;
shrink = [1 1 1];
smooth = [0.1 20];

[Internal, header] = LoadAnalyze('Middle_levelset_Result.hdr', 'Real');
mcubes2VTKFile(Internal, header, vtkfile, threshold, shrink, smooth);


datacolor = [1 1 1];
opacity = 1;
RenderVTKFile(vtkfile, smooth, datacolor, opacity);