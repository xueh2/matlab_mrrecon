
function Rigid_Affine_Hreg_RegistrationRun(target_Dir, target_anatomyDir,...
                                    source_Dir, source_anatomyDir,...
                                    result_Dir, Registration_resultsDir,...
                                    rreg_parameterfile, areg_parameterfile, ...
                                    hreg_parameterfiles, numofparameters)
                                
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

rregname = [nameT '-' nameS '-rreg.dof'];
rregfullname = fullfile(result_Dir, Registration_resultsDir, rregname);

aregname = [nameT '-' nameS '-areg.dof'];
aregfullname = fullfile(result_Dir, Registration_resultsDir, aregname);

hregnames = cell(numofparameters);
for i = 1:numofparameters
    hregnames{i} = [nameT '-' nameS '-hreg' '-' num2str(i) '.dof'];
    hregnames{i} = fullfile(result_Dir, Registration_resultsDir, hregnames{i});
end

% perform the registration

    % rreg
      RigidRegistrationRun2(target, source, rregfullname, rreg_parameterfile);
    
    % areg
      AffineRegistrationRun2(target, source, aregfullname, areg_parameterfile, rregfullname);
    
    % hreg
%     NonRigidRegistrationRun2(target, source, hregnames, hreg_parameterfiles, numofparameters, aregfullname);
return;
