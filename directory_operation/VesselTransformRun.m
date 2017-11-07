
function VesselTransformRun(home, vesselDir, offsetDir, newImageDir, oriImageDir, dofDir, vesselTransformedDir, overlapflag)
% transform vessels into the coordinate of anatomical images
% the vessel extracted is stored as mat in vesselDir, named "vasculature"
% offsetDir stores a txt file recording the leftup and rightdown needed to
% be offset
% dofDir store one dof file transforming vessels into target coordinate

% dof file
filename = fullfile(home, dofDir, '*.dof');
indir = dir(filename);
num_indir = length(indir);
if ( num_indir == 0 )
    disp('can not find dof file ...');
    return;
end

for i = 1:num_indir
    p = strfind(indir(i).name, 'rreg');
    if ( isempty(p) == 0 )
        break;
    end
end

doffilename = fullfile(home, dofDir, indir(i).name);

% load vasculature

% home, extractionResults, tofResampled, tofDir, registrationresults_t2_tof, vesselTransformed_in_t2

filename = fullfile(home, vesselDir, '*.mat');
indir = dir(filename);
num_indir = length(indir);

if ( num_indir == 0 )
    disp('can not find extracted vessels...');
    return
end

vasculature = [];
for k = 1:num_indir
    
    filename = fullfile(home, vesselDir, indir(k).name)
    
    data = load(filename);
    % Validate the MAT-file
    % The file is valid if the variable is called "vasculature" and it has fields called...
    % "count", "vessels" and "totalnumofpoints"
    flds = fieldnames(data);
    num = numel(flds);

    % load vasculature
    index = [];
    for i = 1:num
    % index = find( (flds.name=='vasculature') | ((flds.name=='Vasculature')));
        if ( (strcmp(flds{i}, 'vasculature') == 1) | (strcmp(flds{i}, 'Vasculature') == 1) )
            index = i;
            break;
        end
    end

    if ( isempty(index) == 0 )
        vasculature = eval(['data.' flds{index}])
        break;
    end
end

if ( isempty(vasculature) == 1 )
    disp('no vascualture structure is found...');
    return;
end

% load txt files

offset = [0 0 0 0 0 0];
if ( isempty(offsetDir) == 0 )
  
    filename = fullfile(home, offsetDir, '*.txt');
    indir = dir(filename);
    num_indir = length(indir);
    for i = 1:num_indir
        if ( num_indir == 0 )
            disp('can not find offest file ...');
            break;
        end
    end
    
    name = fullfile(home, offsetDir, indir(1).name);
    
    fid=fopen(name, 'r');
    
    offset = fscanf(fid, '%f');

    fclose(fid);
end

% load new image resolution
filename = fullfile(home, newImageDir, '*.hdr');
indir = dir(filename);
filename = fullfile(home, newImageDir, indir(1).name);
[data, header_new] = LoadAnalyze(filename, 'Grey');

% original image resolution
filename = fullfile(home, oriImageDir, '*.hdr');
indir = dir(filename);
filename = fullfile(home, oriImageDir, indir(1).name);
[data_ori, header_ori] = LoadAnalyze(filename, 'Grey');

% transform the vasculature
offset = offset - 1;
offset_in_Ori_inworld =  [offset(1)*header_ori.xvoxelsize offset(2)*header_ori.yvoxelsize offset(3)*header_ori.zvoxelsize ];

M = eye(4,4);
M(1,1) = header_new.xvoxelsize;% starting from (0 0 0)
M(2,2) = header_new.yvoxelsize;
M(3,3) = header_new.zvoxelsize;
% under the ROI world coordinate
vasculature_in_new_world = vasculatureTransform(vasculature, M)

M = eye(4,4);
M(1,4) = offset_in_Ori_inworld(1);% starting from (0 0 0)
M(2,4) = offset_in_Ori_inworld(2);
M(3,4) = offset_in_Ori_inworld(3);
% under the original world coordinate
vasculature_in_ori_world = vasculatureTransform(vasculature_in_new_world, M)

% ===========================================================
% create the overlap images
if ( overlapflag )
    M = eye(4,4);
    M(1,1) = 1/header_ori.xvoxelsize;% starting from (0 0 0)
    M(2,2) = 1/header_ori.yvoxelsize;
    M(3,3) = 1/header_ori.zvoxelsize;
    % under the original world coordinate
    vasculature_inImage = vasculatureTransform(vasculature_in_ori_world, M);

    volumeSeeds = fuseSeeds_vasculature(header_ori, vasculature_inImage);
    %[overlayMIPxy, overlayMIPyz, overlayMIPzx, overlayVolume] = OverlayRidgePointsOnMIPs(data, volumeSeeds, header_ori, [1 0 0]);
    [overlayMIPxy, overlayMIPyz, overlayMIPzx, overlayVolume] = OverlayRidgePointsOnMIPs(data_ori, volumeSeeds, header_ori, [0 0 1]);
    name = ['overlayMIPxy' '.tiff'];
    filename = fullfile(home, oriImageDir, name);
    imwrite(overlayMIPxy, filename, 'tiff');

    name = ['overlayMIPyz' '.tiff'];
    filename = fullfile(home, oriImageDir, name);
    imwrite(overlayMIPyz, filename, 'tiff');

    name = ['overlayMIPzx' '.tiff'];
    filename = fullfile(home, oriImageDir, name);
    imwrite(overlayMIPzx, filename, 'tiff');
end
% ===========================================================
invertedFlag = 0;
%vasculaturetransformed_in_world = vasculatureDOFtransform(header_ori, vasculature_in_ori_world, doffilename, invertedFlag);
vasculaturetransformed_in_world = vasculatureDOFtransform_CenterDITK(header_ori, vasculature_in_ori_world, doffilename, invertedFlag);

% save the results
filename = fullfile(home, vesselTransformedDir, [vesselTransformedDir '.mat']);
vasculature = vasculaturetransformed_in_world;
save(filename, 'vasculature');

return;
