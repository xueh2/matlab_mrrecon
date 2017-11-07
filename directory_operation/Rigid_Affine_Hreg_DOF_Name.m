
function [rregName, aregName, hregNames] = Rigid_Affine_Hreg_DOF_Name(target_Dir, target_anatomyDir, ...
    source_Dir, source_anatomyDir, numofparameters)
                                
% perform the rigid, affine and non-rigid registration

% targetfile

path_abo = fullfile(target_Dir, target_anatomyDir, '*.hdr' );
indir = dir(path_abo) ;
num = length(indir);
if ( num == 0 )
    disp('empty directory');
    return;
end

target = fullfile(target_Dir, target_anatomyDir, indir(1).name);
[pathstrT,nameT,extT,versnT] = fileparts(target);

% sourcefile

path_abo = fullfile(source_Dir, source_anatomyDir, '*.hdr' );
indir = dir(path_abo);
num = length(indir);
if ( num == 0 )
    disp('empty directory');
    return;
end

source = fullfile(source_Dir, source_anatomyDir, indir(1).name);
[pathstrS,nameS,extS,versnT] = fileparts(source);

% result filenames

rregName = [nameT '-' nameS '-rreg.dof'];

aregName = [nameT '-' nameS '-areg.dof'];

hregNames = cell(numofparameters);
for i = 1:numofparameters
    hregNames{i} = [nameT '-' nameS '-hreg' '-' num2str(i) '.dof'];
end
  
return;
