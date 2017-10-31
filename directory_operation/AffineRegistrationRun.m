
function AffineRegistrationRun(home, subdirectory_target, subdirectory_source, outputDir, parameterfile, dofin)
% rigidly register the image from t2 to tof

path_abo = fullfile(home, subdirectory_target, '*.hdr' );
indir = dir(path_abo) ;
num = length(indir);
if ( num == 0 )
    disp('empty directory');
    return;
end

path_abo_source = fullfile(home, subdirectory_source, '*.hdr' );
indir_source = dir(path_abo_source) ;
num_source = length(indir_source);
if ( num_source == 0 )
    disp('empty directory');
    return;
end

% only register the first found target file and source file with areg
target = fullfile(home, subdirectory_target, indir(1).name);
[pathstrT,nameT,extT,versnT] = fileparts(target);
source = fullfile(home, subdirectory_source, indir_source(1).name);
[pathstrS,nameS,extS,versnT] = fileparts(source);

dofname = [nameT '_' nameS '_areg' '.dof'];
dofname = fullfile(home, outputDir, dofname);

command = ['areg' ' ' target ' ' source ' '  '-dofout' ' ' dofname];
if ( isempty(parameterfile) == 0 )
    command = [command ' '  '-parameter' ' ' parameterfile];
end

if ( isempty(dofin) == 0 )
    command = [command ' ' '-dofin' ' ' dofin];
end
[s, w] = dos(command, '-echo');
return