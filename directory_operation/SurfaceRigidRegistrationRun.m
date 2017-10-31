
function SurfaceRigidRegistrationRun(target, source, dofname, locator, iterations, invert)
% <-locator>          Locator: 0 = cell locator, 1 = point locator, 2 = kd-tree locator (default = 1) 
% <-dofin name>       Name of input file 
% <-dofout name>      Name of output file 
% <-epsilon>          Value for espilon (default=0.01) 
% <-clean>            Clean polydata (default OFF) 
% <-symmetric>        Use symmetric distance (default OFF) 
% <-ignoreedges>      Ignores edges in ICP (default OFF) 
% <-invert>           Save the inverse transformation (default OFF) 
% <-iterations>       Number of 3D registration iterations (default 100) 

command = ['srreg' ' ' target ' ' source ' ' '-dofout' ' ' dofname];

if ( isempty(locator) == 0 )
    command = [command ' ' '-locator' ' ' num2str(locator)];
end

if ( isempty(iterations) == 0 )
    command = [command ' ' '-iterations' ' ' num2str(iterations)];
end

if ( invert == 1 )
    command = [command ' ' '-invert'];
end

disp(command)
[s, w] = dos(command, '-echo');
return