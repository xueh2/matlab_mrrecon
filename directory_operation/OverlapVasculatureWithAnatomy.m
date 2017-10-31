
function OverlapVasculatureWithAnatomy(home, vesselTransformedDir, anatomyDir, default_radius)

% the vessel has been transformed into the world coordinate of anatomy, following
% the DITK coordinate convention
% vessel should be mapped into the
% leftup-corner coordinate from its original position in
% vesselTransformedDir

% overlap the vasculature on the anatomy as white dots

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

% plotTwoVasculature(vasculature_updown, vasculature_in_leftup_world, 'b', 'r')

% to the image coordinate
M = eye(4,4);
M(1,1) = 1/header.xvoxelsize;% starting from (0 0 0)
M(2,2) = 1/header.yvoxelsize;
M(3,3) = 1/header.zvoxelsize;
vasculature_in_leftup_image = vasculatureTransform(vasculature_in_leftup_world, M)

% plotTwoVasculature(vasculature_in_leftup_world, vasculature_in_leftup_image, 'b', 'r')

% mapping
default_radius = (header.xvoxelsize + header.yvoxelsize + header.zvoxelsize ) / 3;
disp('mapping ... ');
volume = fuseSeeds_vasculature_radius(header, vasculature_in_leftup_image, default_radius);

maxI = max(max(max(data)));
data( find(volume == true) ) = 2 * maxI;

% save into vesselTransformedDir
filename = fullfile(home, anatomyDir, '*.hdr');
indir = dir(filename);
filename = ['vessel_' indir(1).name];

filename = fullfile(home, vesselTransformedDir, filename);
SaveAnalyze(data, header, filename, 'Grey');

disp('mapping finished ... ');

return

