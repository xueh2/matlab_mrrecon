
function [data, g] = CreateDigitalSphere(header, radius)
% create a digial sphere
% the resolution should be isotropic and the radius is in the unit of mm

xvoxelsize = header.xvoxelsize;
yvoxelsize = header.yvoxelsize;
zvoxelsize = header.zvoxelsize;
xsize = header.xsize;
ysize = header.ysize;
zsize = header.zsize;

%--------------------------------------------------------------------------
% Create the grid.
g = [];
g.dim = 3;
% ==========================
g.min = [-yvoxelsize * (ysize-1)/2.0; -xvoxelsize * (xsize-1)/2.0; -zvoxelsize * (zsize-1)/2.0];
g.max = [yvoxelsize * (ysize-1)/2.0; xvoxelsize * (xsize-1)/2.0; zvoxelsize * (zsize-1)/2.0];

g.dx = [yvoxelsize; xvoxelsize; zvoxelsize];

g.bdry = @addGhostExtrapolate;

g = processGrid(g);

%--------------------------------------------------------------------------

data = shapeSphere(g, [0 0 0], radius);

return;