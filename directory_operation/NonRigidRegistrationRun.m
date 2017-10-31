
function NonRigidRegistrationRun(home_target, subdirectory_target,...
                                    home_source, subdirectory_source,...
                                    home_output, outputDir,...
                                    parameterfile, numofparameters, dofin)
% rigidly register the image from t2 to tof

path_abo = fullfile(home_target, subdirectory_target, '*.hdr' );
indir = dir(path_abo) ;
num = length(indir);
if ( num == 0 )
    disp('empty directory');
    return;
end

path_abo_source = fullfile(home_source, subdirectory_source, '*.hdr' );
indir_source = dir(path_abo_source) ;
num_source = length(indir_source);
if ( num_source == 0 )
    disp('empty directory');
    return;
end

% only register the first found target file and source file with areg
target = fullfile(home_target, subdirectory_target, indir(1).name);
[pathstrT,nameT,extT,versnT] = fileparts(target);
source = fullfile(home_source, subdirectory_source, indir_source(1).name);
[pathstrS,nameS,extS,versnT] = fileparts(source);

command = ['hreg' ' ' target ' ' source ];
if ( isempty(parameterfile) == 0 )
    
    command = [command ' '  '-parameter' num2str(numofparameters)];
    
    for i = 1:numofparameters
        command = [command ' ' parameterfile{i}];
    end
end

if ( isempty(dofin) == 0 )
    command = [command ' ' '-dofin' ' ' dofin];
end

command = [command ' '  '-dofout' ];
for i = 1:numofparameters
    dofname = [nameT '_' nameS '_hreg' '_' num2str(i) '.dof'];
    dofname = fullfile(home_output, outputDir, dofname);
    command = [command ' ' dofname];
end
if ( numofparameters == 0 )
    dofname = [nameT '_' nameS '_hreg' '.dof'];
    dofname = fullfile(home_output, outputDir, dofname);
    command = [command ' ' dofname];
end

[s, w] = dos(command, '-echo');
return