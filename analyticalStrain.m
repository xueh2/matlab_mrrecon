function [rad_strain, circ_strain, sheer] = analyticalStrain(mask, de_all, dr_all)

imsz = size(dr_all);
if (numel(imsz) == 2)
    slices = 1;
else
    slices = imsz(3);
end

rad_strain = zeros(imsz);
circ_strain = zeros(imsz);
sheer = zeros(imsz);

[Emesh, Rmesh] = meshgrid(1:imsz(2), 1:imsz(1));
% find centroid of mask
centroidE = sum(sum(Emesh(mask == 1)))/sum(sum(mask));
centroidR = sum(sum(Rmesh(mask == 1)))/sum(sum(mask));
centroid = [centroidE, centroidR];

for i = 1:slices
    dr = dr_all(:,:,i);
    de = de_all(:,:,i);
    
    Y = -(Rmesh - centroidR);
    X = -(centroidE - Emesh);

    [ddrdr, ddrde] = imgradientxy_cine(dr, 'central');
    [ddedr, ddede] = imgradientxy_cine(-de, 'central');

    F00 = (1+ddrdr);
    F01 = ddrde;
    F10 = ddedr;
    F11 = (1+ddede);
    
    % 1/2 [F00 F10] [F00 F01]  - [1 0]
    %     [F01 F11] [F10 F11]    [0 1]

    E00 = 1/2*(F00.^2 + F10.^2 - 1).*mask;
    E01 = 1/2*(F00.*F01 + F10.*F11).*mask;
    E10 = 1/2*(F00.*F01 + F10.*F11).*mask;
    E11 = 1/2*(F01.^2 + F11.^2 - 1).*mask;
 
    thetas = atan(Y./(X+1e-8)) + pi*(X < 0) + pi*2*(X>=0).*(Y <0);
    rad_strain(:, :, i) = E00.*cos(thetas).^2 + E11.*sin(thetas).^2 ...
        + 2*E01.*sin(thetas).*cos(thetas);
    circ_strain(:, :, i) = E00.*sin(thetas).^2 + E11.*cos(thetas).^2 ...
        - 2*E01.*sin(thetas).*cos(thetas);
    sheer(:, :, i) = (E11-E00).*sin(thetas).*cos(thetas) ...
        + E01.*(cos(thetas).^2 - sin(thetas).^2);
    
end
end