function [C, rho, S] = CoilSensitivity_Souheil_Walsh_PowerMethod(data, kSize)
% Modified Walsh method
%  Author: Souheil Inati, NIH
%
% [C, rho, S] = CoilSensitivity_Souheil_Walsh(data, kSize)
%
% Input:
%   data: (Nx,Ny,Nz) array
%     2D complex images (Nx,Ny) from Nc coils
%     assumed to be prewhitened
%
% Output:
%   rho: (Nx,Ny)
%     adaptively combined image
%
%   C: (Nx,Ny,Nc)
%     estimated relative coil sensitivities
%
%   S: (Nx,Ny)
%     smooth estimate of sqrt(Sum_j |Cj|^2 * |rho|^2)
%
%

[Nx,Ny,Nc] = size(data);

% Block size
% These have to be odd for the code below
Bx = kSize(1); By = kSize(2);

% Usefule shortcut parameters
% bxmin = floor(Bx/2)+1; bymin = floor(By/2)+1;
% bxmax = Nx-floor(Bx/2); bymax = Ny-floor(By/2);
dx = floor(Bx/2); dy = floor(By/2);

% Initialize the result
C   = zeros(Nx,Ny,Nc);  % sensitivities (normalized)
S   = zeros(Nx,Ny);     % |C|*|rho|
rho = zeros(Nx,Ny);     % adaptive combination

% Loop over the blocks
% We ignore the edges here for simplicity
% A full implementaion would handle the edges
% by zeros, circular, or mirror symmetry.
% for x = bxmin:bxmax
%     for y = bymin:bymax
%         
%         D = reshape(data(x-dx:x+dx,y-dy:y+dy,:),[Bx*By,Nc]);
for x = 1:Nx
    for y = 1:Ny
        
        indX = x-dx:x+dx;
        indY = y-dy:y+dy;
        
        p = find(indX<1);
        if ( ~isempty(p) )
            indX(p) = Nx - indX(p);
        end
        
        p = find(indY<1);
        if ( ~isempty(p) )
            indY(p) = Ny - indY(p);
        end
        
        p = find(indX>Nx);
        if ( ~isempty(p) )
            indX(p) = -Nx + indX(p);
        end
        
        p = find(indY>Ny);
        if ( ~isempty(p) )
            indY(p) = -Ny + indY(p);
        end        
        
        D = reshape(data(indX,indY,:),[Bx*By,Nc]);
        
        % The SVD way
        % [U,S,V]   = svd(D);
        % u1 = U(:,1); % first left eigenvector (Bx*By,1)
        % v1 = V(:,1); % first right eigenvector (Nc,1)
        % s1 = S(1,1); % first singular value
        
        % The power method way
        % 3 iterations
        % we solve for v1
        % initialize to the mean of the data
        v1 = transpose(mean(D,1)); 
        v1 = v1/norm(v1);
        for iter = 1:3
            v1 = D'*D*v1; 
            v1 = v1/norm(v1);
        end
        % compute u1 and s1
        u1 = D*v1;
        s1 = norm(u1);
        u1 = u1/s1;

        % Compute the phase of the average of u1
        mu1 = mean(u1);
        theta = angle(mu1);
        
        % Smoothed coil sensitivity magnitude * rho magnitude
        % Normalize the singular value by the box size
        S(x,y) = s1/sqrt(Bx*By);
        
        % Coil sensitivies
        % Put the phase into the coil sensitivities
        C(x,y,:)  = exp(1i*theta)*conj(v1);
        
        % Combined image
        % Remove the phase from the combined image
        rho(x,y)  = u1(floor(Bx*By/2)+1)*exp(-1i*theta)*s1;        

    end % Loop over y
end % Loop over x

end