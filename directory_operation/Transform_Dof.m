
function Transform_Dof(source_Dir, source_anatomyDir,...
                       output_Dir, output_anatomyDir, outputfile, ...
                       dof_Dir, dof_anatomyDir)
                                
% perform the rigid registration
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

% outputfile
output = fullfile(output_Dir, output_anatomyDir, outputfile );

% doffile

path_abo = fullfile(dof_Dir, dof_anatomyDir, '*.dof' );
indir = dir(path_abo);
num = length(indir);
if ( num == 0 )
    disp('empty directory');
    return;
end
dof = fullfile(dof_Dir, dof_anatomyDir, indir(1).name);

% transformation
TransformationRun(source, output, dof);
    
return;
