
function [shockpoints, shockvalues, shocks] = GetShockPoints(grid, TTR, SignedPressureForce_GMEnhanced, SDF, shockThreshold)
% get the shock points
% shockpoints: N*1 vector, indexes of shock points
% shockvalues: correponding shock values F(x)*norm(gradient(D(x)))

deriv_X = centeredFirstSecond(grid, TTR, 2);
deriv_y = centeredFirstSecond(grid, TTR, 1);
deriv_z = centeredFirstSecond(grid, TTR, 3);

grad_Mag = sqrt(deriv_X.^2+deriv_y.^2+deriv_z.^2);
clear deriv_X deriv_Z deriv_Y

shocks = SignedPressureForce_GMEnhanced .* grad_Mag;

shockpoints = find( (SDF>0) & (shocks<=shockThreshold) );
shockvalues = shocks(shockpoints);

return;