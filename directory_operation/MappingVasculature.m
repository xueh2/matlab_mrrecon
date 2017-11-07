
function MappingVasculature(targetDir, mappingresultsDir, interReg_Dir, target_anatomyDir, ...
    sourceDir, vessel_in_worldDir, source_anatomyDir)

% mapping vasculature in sourceDir/vessel_in_worldDir into targetDir/mappingresultsDir.
% dof file is stored in targetDir/interReg_Dir
% dof file name: nameT-nameS-rreg/areg/hreg-num.dof
% mapping results: vasculature-nameS-rreg/areg/hreg-num.mat

% load vasculature
filename = fullfile(sourceDir, vessel_in_worldDir, '*.mat');
indir = dir(filename);
num_indir = length(indir);

if ( num_indir == 0 )
    disp('can not find extracted vessels (mat file)...');
    return
end

vasculature = [];
for k = 1:num_indir
    filename = fullfile(sourceDir, vessel_in_worldDir, indir(k).name);
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
        vasculature_in_ori_world = eval(['data.' flds{index}]);
        break;
    end
end

if ( isempty(vasculature_in_ori_world) == 1 )
    disp('no vascualture structure is found...');
    return;
end

% dof file

filename = fullfile(targetDir, interReg_Dir, '*.dof');
indir = dir(filename);
num_indir = length(indir);
if ( num_indir == 0 )
    disp('can not find extracted vessels (mat file)...');
    return
end

place = find(targetDir(length(targetDir):-1:1) == filesep );
ctarget = targetDir(length(targetDir)-place(1)+2:length(targetDir));

place = find(sourceDir(length(sourceDir):-1:1) == filesep );
csource = sourceDir(length(sourceDir)-place(1)+2:length(sourceDir));

doffiles = cell(0);
sources = cell(0);
rahs = cell(0);
numofhs = cell(0);
for i = 1:num_indir
    
    [target, source, rah, numofh] = ParseDofName(indir(i).name);
    
    if ( (strcmp(ctarget, target) == 1) & (strcmp(csource, source) == 1) )
        doffiles = [doffiles {indir(i).name}];
        sources = [sources {source}];
        rahs = [rahs {rah}];
        numofhs = [numofhs {numofh}];
    end
end

% header_target
filename = fullfile(targetDir, target_anatomyDir, '*.hdr');
indir = dir(filename);
num_indir = length(indir);
if ( num_indir == 0 )
    disp('can not find anatomy (hdr file)...');
    return
end

filename = fullfile(targetDir, target_anatomyDir, indir(1).name);
[data, header_target] = LoadAnalyze(filename, 'Grey');

% header_source
filename = fullfile(sourceDir, source_anatomyDir, '*.hdr');
indir = dir(filename);
num_indir = length(indir);
if ( num_indir == 0 )
    disp('can not find anatomy (hdr file)...');
    return
end

filename = fullfile(sourceDir, source_anatomyDir, indir(1).name);
[data, header_source] = LoadAnalyze(filename, 'Grey');


% transform vessels
[lines, num] = size(doffiles);
for i = 1:num
    
    source = sources{i};
    
    if ( strcmp(rahs{i}, 'rreg') == 1 )
        
        name = ['vasculature-' source '-' rahs{i} '.mat'];
        fullname = fullfile(targetDir, mappingresultsDir, name);
        doffullname = fullfile(targetDir, interReg_Dir, doffiles{i});
        
        invertedFlag = 0;
        %vasculatureTransformed = vasculatureDOFtransform_twoHeaders(header_source, header_target, vasculature_in_ori_world, doffullname, invertedFlag);
%         vasculatureTransformed = vasculatureDOFtransform_twoHeaders(header_source, header_source,...
%             vasculature_in_ori_world, doffullname, invertedFlag);
        
        vasculatureTransformed = vasculatureDOFtransform_Already_CenterDITK(vasculature_in_ori_world, doffullname, invertedFlag);
        
    end
    
    if ( strcmp(rahs{i}, 'areg') == 1 )
        
        name = ['vasculature-' source '-' rahs{i} '.mat'];
        fullname = fullfile(targetDir, mappingresultsDir, name);
        doffullname = fullfile(targetDir, interReg_Dir, doffiles{i});
        
        invertedFlag = 0;
        %vasculatureTransformed = vasculatureDOFtransform_twoHeaders(header_source, header_target, vasculature_in_ori_world, doffullname, invertedFlag);
%         vasculatureTransformed = vasculatureDOFtransform_twoHeaders(header_source, header_source,...
%             vasculature_in_ori_world, doffullname, invertedFlag);
        
        vasculatureTransformed = vasculatureDOFtransform_Already_CenterDITK(vasculature_in_ori_world, doffullname, invertedFlag);

    end
    
    if ( strcmp(rahs{i}, 'hreg') == 1 )
        
        name = ['vasculature-' source '-' rahs{i} '-' num2str(numofhs{i}) '.mat'];
        fullname = fullfile(targetDir, mappingresultsDir, name);
        doffullname = fullfile(targetDir, interReg_Dir, doffiles{i});
        
        invertedFlag = 0;
        %vasculatureTransformed = vasculatureDOFtransform_twoHeaders(header_source, header_target, vasculature_in_ori_world, doffullname, invertedFlag);
%         vasculatureTransformed = vasculatureDOFtransform_twoHeaders(header_source, header_source,...
%             vasculature_in_ori_world, doffullname, invertedFlag);

        vasculatureTransformed = vasculatureDOFtransform_Already_CenterDITK(vasculature_in_ori_world, doffullname, invertedFlag);
    end
    
    disp(fullname);
    disp('==============================================');
    vasculature = vasculatureTransformed;
    save(fullname, 'vasculature');
    
end

return;