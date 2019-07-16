function [rad_strain, circ_strain, sheer] = analyticalStrain(mask, dr_all, de_all)

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
    
    X = Rmesh - centroidR;
    Y = centroidE - Emesh;

    [ddrde, ddrdr] = imgradientxy_cine(dr, 'central');
    [ddede, ddedr] = imgradientxy_cine(-de, 'central');

    F00 = 1 + ddrdr;
    F01 = ddrde;
    F10 = ddedr;
    F11 = 1 + ddede;
    
    
    
%     figure
%     imagescn([F00.*mask, F01.*mask; F10.*mask, F11.*mask])
    
    % 1/2 [F00 F10] [F00 F01]  - [1 0]
    %     [F01 F11] [F10 F11]    [0 1]

    E00 = 1/2*F00.^2 + 1/2*F10.^2 - 1;
    E01 = 1/2*F00.*F01 + 1/2*F10.*F11;
    E10 = 1/2*F00.*F01 + 1/2*F10.*F11;
    E11 = 1/2*F01.^2 + 1/2*F11.^2 - 1;

%     figure
%     imagescn([E00.*mask, E01.*mask; E10.*mask, E11.*mask])
    
    thetas = atan(Y./X) + pi*(X < 0) + pi*2*(X>=0).*(Y <0);
%     figure
%     imagescn(thetas)
    sheer(:, :, i) = 2*(E11-E00).*sin(thetas).*cos(thetas) ...
        + 2*E01.*(cos(thetas).^2 - sin(thetas).^2);
    rad_strain(:, :, i) = E00.*cos(thetas).^2 + E11.*sin(thetas).^2 ...
        + 2*E01.*sin(thetas).*cos(thetas);
    thetas = thetas + pi/2;
    circ_strain(:, :, i) = E00.*cos(thetas).^2 + E11.*sin(thetas).^2 ...
        + 2*E01.*sin(thetas).*cos(thetas);

end
end