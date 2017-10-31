
function RigidRegistrationRun(home, subdirectory_tof, subdirectory_anatomy, outputDir, parameterfile)
% rigidly register the image from t2 to tof

path_abo = fullfile(home, subdirectory_tof, '*.hdr' );
indir = dir(path_abo) ;
num = length(indir);
if ( num == 0 )
    disp('empty directory');
    return;
end

path_abo_t2 = fullfile(home, subdirectory_anatomy, '*.hdr' );
indir_t2 = dir(path_abo_t2) ;
num_t2 = length(indir_t2);
if ( num_t2 == 0 )
    disp('empty directory');
    return;
end

% only register the first found tof file and t2 file with rreg
target = fullfile(home, subdirectory_tof, indir(1).name);
[pathstrT,nameT,extT,versnT] = fileparts(target);
source = fullfile(home, subdirectory_anatomy, indir_t2(1).name);
[pathstrS,nameS,extS,versnT] = fileparts(source);

dofname = [nameT '_' nameS '_rreg' '.dof'];
dofname = fullfile(home, outputDir, dofname);
if ( isempty(parameterfile) == 0 )
    command = ['rreg' ' ' target ' ' source ' ' '-parameter' ' ' parameterfile ' '  '-dofout' ' ' dofname];
else
    command = ['rreg' ' ' target ' ' source ' ' '-dofout' ' ' dofname];
end
[s, w] = dos(command, '-echo');
return