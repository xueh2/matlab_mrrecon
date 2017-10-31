
% home = 'D:\work\directory_operation';
% subdirectory = 'sense';
% suffix = 'tof';
% 
% % FormateFileNames(home, subdirectory, suffix);
% % 
% %  [subdir, num] = FindAllDirectory(home)
% %  
%  resolution = [0.2 0.2 0.2];
%  
%  ResampleRun(home, subdirectory, resolution)
 
home = 'D:\work\hui';
subdirectory_tof = 'sense';
subdirectory_anatomy = 't2';

CreateDirs(home, 'test')

ResampleRun(home, subdirectory_anatomy, [0.2 0.2 0.2])

% FormateFileNames(home, subdirectory_tof, 'tof');
% FormateFileNames(home, subdirectory_anatomy, 't2');
%  
%  RigidRegistrationRun(home, subdirectory_tof, subdirectory_anatomy, [])
 
 target = 'D:\work\hui\sense\hui_tof_1.hdr';
 source = 'D:\work\hui\t2\hui_t2_1.hdr';
 dofin = 'D:\work\hui\hui_tof_1_hui_t2_1_areg.dof';
 rview(target, source, dofin)
 
home = 'D:\work\hui';
subdirectory_tof = 'sense';
subdirectory_anatomy = 't2';
dofin = 'D:\work\hui\hui_tof_1_hui_t2_1_rreg.dof';
 AffineRegistrationRun(home, subdirectory_tof, subdirectory_anatomy, [], dofin)
 
 dofin = 'D:\work\hui\hui_tof_1_hui_t2_1_areg.dof';
 NonRigidRegistrationRun(home, subdirectory_tof, subdirectory_anatomy, [], 0, dofin)