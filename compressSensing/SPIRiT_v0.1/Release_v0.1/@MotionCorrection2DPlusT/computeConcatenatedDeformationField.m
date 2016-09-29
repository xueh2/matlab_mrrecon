function a = computeConcatenatedDeformationField(a, dx, dy, dxInv, dyInv)
% compute the concatenated deformation field
% x_new = x + dx + a.dx(dx+x, dy+y)
% y_new = y + dy + a.dy(dx+x, dy+y)
% xInv_new = x + a.dxInv + dxInv(a.dxInv+x, a.dyInv+y)
% yInv_new = y + a.dyInv + dyInv(a.dxInv+x, a.dyInv+y)

if isa(a,'MotionCorrection2DPlusT') == 0
    error('In computeConcatenatedDeformationField(a, dx, dy, dxInv, dyInv), a must be MotionCorrection2DPlusT operator');
end

[Nfe, Npe, Nrep] = size(dx);

dxNew = dx;
dyNew = dy;
dxInvNew = dxInv;
dyInvNew = dyInv;

[X, Y] = meshgrid(0:Npe-1, 0:Nfe-1);

for r=1:Nrep
    
    % compute the dxNew, dyNew
    xTemp = dx(:,:,r) + X;
    yTemp = dy(:,:,r) + Y;
    
    dxTemp = interp2(X, Y, a.dx(:,:,r), xTemp, yTemp, 'linear');
    dyTemp = interp2(X, Y, a.dy(:,:,r), xTemp, yTemp, 'linear');
    
    dxNew(:,:,r) = dx(:,:,r) + dxTemp;
    dyNew(:,:,r) = dy(:,:,r) + dyTemp;
    
    % compute the dxInvNew, dyInvNew
    xInvTemp = a.dxInv(:,:,r) + X;
    yInvTemp = a.dyInv(:,:,r) + Y;

    dxInvTemp = interp2(X, Y, dxInv(:,:,r), xInvTemp, yInvTemp, 'linear');
    dyInvTemp = interp2(X, Y, dyInv(:,:,r), xInvTemp, yInvTemp, 'linear');
    
    dxInvNew(:,:,r) = a.dxInv(:,:,r) + dxInvTemp;
    dyInvNew(:,:,r) = a.dyInv(:,:,r) + dyInvTemp;
end

a.dx = single(dxNew);
a.dy = single(dyNew);
a.dxInv = single(dxInvNew);
a.dyInv = single(dyInvNew);
