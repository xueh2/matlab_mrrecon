
function RenderVasculatureWithAnatomy(home, vesselTransformedDir, anatomyDir)

% the vessel has been transformed into the world coordinate of anatomy, following
% the DITK coordinate convention
% to render them correctly using VTK, vessel should be transformed into the
% leftup-corner coordinate from its original position in
% vesselTransformedDir

% load vasculature
filename = fullfile(home, vesselTransformedDir, '*.mat');
indir = dir(filename);
num_indir = length(indir);

if ( num_indir == 0 )
    disp('can not find transformed vessels...');
    return
end

vasculature = [];
for k = 1:num_indir
    
    filename = fullfile(home, vesselTransformedDir, indir(k).name)
    
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

% anatomy images
filename = fullfile(home, anatomyDir, '*.hdr');
indir = dir(filename);
filename = fullfile(home, anatomyDir, indir(1).name);
[data, header] = LoadAnalyze(filename, 'Grey');

% transform the vasculature into the leftup-corner coordinate

M = eye(4,4);
M(1,1) = 1;% starting from (0 0 0)
M(2,2) = -1;
M(3,3) = 1;
vasculature_updown = vasculatureTransform(vasculature, M);

% plotTwoVasculature(vasculature, vasculature_updown, 'b', 'r')

M = eye(4,4);
M(1,4) = header.xvoxelsize * (header.xsize-1)/2.0;% starting from (0 0 0)
M(2,4) = header.yvoxelsize * (header.ysize-1)/2.0;
M(3,4) = header.zvoxelsize * (header.zsize-1)/2.0;
vasculature_in_leftup_world = vasculatureTransform(vasculature_updown, M);

M = eye(4,4);
M(1,1) = 1/header.xvoxelsize;% starting from (0 0 0)
M(2,2) = 1/header.yvoxelsize;
M(3,3) = 1/header.zvoxelsize;
vasculature_in_leftup_image = vasculatureTransform(vasculature_updown, M);

% plotTwoVasculature(vasculature_updown, vasculature_in_leftup_world, 'b', 'r')

% to the image coordinate
% M = eye(4,4);
% M(1,1) = 1/header.xvoxelsize;% starting from (0 0 0)
% M(2,2) = 1/header.yvoxelsize;
% M(3,2) = 1/header.zvoxelsize;
% vasculature_in_leftup_world = vasculatureTransform(vasculature_updown, M)

% rendering...
rcolor = createColor(vasculature_in_leftup_world, [0.8 0 0]);
opacity = createOpacity(vasculature_in_leftup_world, 0.75);
if ( isempty(vasculature.vessels{1}.radius) == 1 )
    RenderVasculature(vasculature_in_leftup_world, 0.05, rcolor, opacity, data, header);
else
    RenderVasculature(vasculature_in_leftup_world, -1, rcolor, opacity, data, header);
end

return

