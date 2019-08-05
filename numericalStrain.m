function [radial_strain, circumfrential_strain] = numericalStrain(mask, dx_all, dy_all, step_size)

if(nargin<4)
    step_size = 0.5;
end

imsz = size(dx_all);
if (numel(imsz) == 2)
    slices = 1;
else
    slices = imsz(3);
end

if (numel(size(mask)) == 2)
    mask_slices = 1;
else
    mask_slices = imsz(3);
end

radial_strain = zeros(imsz);
circumfrential_strain = zeros(imsz);

centroidE = zeros(mask_slices, 1);
centroidR = zeros(mask_slices, 1);
for i = 1:mask_slices
    [Emesh, Rmesh] = meshgrid(1:imsz(2), 1:imsz(1));
    % find centroid of mask
    if mask_slices == 1
        centroidE(i) = sum(sum(Emesh(mask == 1)))/sum(sum(mask));
        centroidR(i) = sum(sum(Rmesh(mask == 1)))/sum(sum(mask));
    else
        centroidE(i) = sum(sum(Emesh(mask(:, :, i) == 1)))/sum(sum(mask(:, :, i)));
        centroidR(i) = sum(sum(Rmesh(mask(:, :, i) == 1)))/sum(sum(mask(:, :, i)));
    end
    centroid = [centroidR(i), centroidE(i)];
end
    
    
for i = 1:slices
    dx = dx_all(:,:,i);
    dy = dy_all(:,:,i);
    if mask_slices == 1
        X = Rmesh - centroidR;
        Y = centroidE - Emesh;
    else
        X = Rmesh - centroidR(i);
        Y = centroidE(i) - Emesh;
    end

    % get radial thetas
    theta = atan(Y./(X+1e-8)) + pi*(X < 0) + pi*2*(X>=0).*(Y <0);

    % get endpoints
    X_in = Rmesh + step_size*cos(theta);
    X_out = Rmesh - step_size*cos(theta);
    Y_in = Emesh - step_size*sin(theta);
    Y_out = Emesh + step_size*sin(theta);

    % get distances
    distances = 2*step_size*ones(imsz(1), imsz(2));

    % get interpolated points
    dx_in = interp2(dx, Y_in, X_in);
    dx_out = interp2(dx, Y_out, X_out);
    dy_in = interp2(dy, Y_in, X_in);
    dy_out = interp2(dy, Y_out, X_out);

    X_prime_in = dx_in + X_in;
    X_prime_out = dx_out + X_out;
    Y_prime_in = dy_in + Y_in;
    Y_prime_out = dy_out + Y_out;

    % get radial projection magnitude
    a_x = X_prime_out - X_prime_in;
    a_y = Y_prime_out - Y_prime_in;
    b_x = X_out - X_in;
    b_y = Y_out - Y_in;
    comp_ab_rad = a_x.*b_x + a_y.*b_y;

    % calculate radial_strain
    radial_strain(:, :, i) = (comp_ab_rad - distances)./distances;

    % get rotated thetas
    theta_rot = theta + pi/2;

    % get circ endpoints
    X_in = Rmesh + step_size*cos(theta_rot);
    X_out = Rmesh - step_size*cos(theta_rot);
    Y_in = Emesh - step_size*sin(theta_rot);
    Y_out = Emesh + step_size*sin(theta_rot);

    % get distances
%     distances = 2*step_size* ones(imsz(1), imsz(2));

    % get circ interpolated points
    dx_in = interp2(dx, Y_in, X_in);
    dx_out = interp2(dx, Y_out, X_out);
    dy_in = interp2(dy, Y_in, X_in);
    dy_out = interp2(dy, Y_out, X_out);

    X_prime_in = dx_in + X_in;
    X_prime_out = dx_out + X_out;
    Y_prime_in = dy_in + Y_in;
    Y_prime_out = dy_out + Y_out;

    % get circ projection magnitude
    a_x_rot = X_prime_out - X_prime_in;
    a_y_rot = Y_prime_out - Y_prime_in;
    b_x_rot = X_out - X_in;
    b_y_rot = Y_out - Y_in;
    comp_ab_circ = a_x_rot.*b_x_rot + a_y_rot.*b_y_rot;

    % calculate circumfrential strain
    circumfrential_strain(:, :, i) = (comp_ab_circ - distances)./distances;
end

ind = find(isnan(radial_strain));
radial_strain(ind) = 0;

end