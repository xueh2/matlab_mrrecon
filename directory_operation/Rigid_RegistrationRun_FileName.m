
function Rigid_RegistrationRun_FileName(target_FullName, source_FullName, dofout_FullName, rreg_parameterfile)
                                
% perform the rigid registration

% targetfile
target = target_FullName;
[pathstrT,nameT,extT,versnT] = fileparts(target);

% sourcefile
source = source_FullName;
[pathstrS,nameS,extS,versnT] = fileparts(source);

% result filenames
rregfullname = dofout_FullName;

% perform the registration

% rreg
RigidRegistrationRun2(target, source, rregfullname, rreg_parameterfile);
    
return;
