
% parameters

% NumofLevels <-- the level nums of iterations
% maximal_inflation_factor [0 1]
% NumberOfIterations = 5; <-- the inflation iteration between two checking
% RelaxationFactor = 0.5; < lemda in inflation
% ThresRatio = 0.01; <-- the threshold to stop
% numofSurfaceNonRigid

% targetImageName
% Internal_Target_Name
% 
% sourceImageName
% Internal_Source_Name

% CortexRegistationDir
cd(CortexRegistationDir);

locator = 2;
iterations = 250;
epsilon = 0.005;
symmetric = [];

numofparameters = 5; % 20mm --> 1.25mm

spacing = 20;

if ( exist('Prefix') == 0 )
    Prefix = '';
end

if ( exist('numofSurfaceNonRigid') == 0 )
    numofSurfaceNonRigid = 4;
end

% -------------------------------------------------------------------------
surfacefile_Target = Internal_Target_Name;
vtkfile_Target = 'internal_target.vtk';
vtkfile_Target_inflated = 'internal_target_inflated.vtk';
image_Target = targetImageName;

surfacefile_Source = Internal_Source_Name;
vtkfile_Source = 'internal_source.vtk';
vtkfile_Source_inflated = 'internal_source_inflated.vtk';
image_Source = sourceImageName;


threshold = 0;
shrink = [1 1 1];
smooth = [0.1 20];

datacolor = [1 1 1];
opacity = 1;

rreg_name = 'image_rigid.dof';
areg_name = 'image_affine.dof';

% -------------------------------------------------------------------------
% image registration
target = image_Target;
source = image_Source;

% all subject directories
rreg_parameterfile = 'K:\more_neonates_images_LS\parameters\parameters.rreg';
areg_parameterfile = 'K:\more_neonates_images_LS\parameters\parameters.areg';

numofparameters = 4; % 20mm --> 1.25mm
hreg_parameterfiles = cell(numofparameters);
hreg_parameterfiles{1} = 'K:\more_neonates_images_LS\parameters\parameters-20mm.mreg';
hreg_parameterfiles{2} = 'K:\more_neonates_images_LS\parameters\parameters-10mm.mreg';
hreg_parameterfiles{3} = 'K:\more_neonates_images_LS\parameters\parameters-5mm.mreg';
hreg_parameterfiles{4} = 'K:\more_neonates_images_LS\parameters\parameters-2.5mm.mreg';
hreg_parameterfiles{5} = 'K:\more_neonates_images_LS\parameters\parameters-2.5mm.mreg';

dofnames = cell(numofparameters);
dofnames{1} = 'image_nonrigid_20mm.dof';
dofnames{2} = 'image_nonrigid_10mm.dof';
dofnames{3} = 'image_nonrigid_5mm.dof';
dofnames{4} = 'image_nonrigid_2.5mm.dof';
dofnames{5} = 'image_nonrigid_1.25mm.dof';

controlPoints = [20 20 20];
if ( isempty(dir(rreg_name))==1 )
    RigidRegistrationRun2(target, source, rreg_name, rreg_parameterfile);
end

if ( isempty(dir(areg_name))==1 )
    AffineRegistrationRun2(target, source, areg_name, areg_parameterfile, rreg_name);
end
NonRigidRegistrationRun2(target, source, dofnames, hreg_parameterfiles, numofparameters, areg_name, 0, controlPoints);
%------------------------------------------------------------------------