
function csm = CoilSensitivity_Souheil_Walsh(data, kSize, plotFlag)

header = CreateFtkHeaderInfo(data, [1 1 1 1]);
csm = Matlab_PerformCoilMapEstimation( single(data), header, 'Souheil', 3, kSize(1), 90, size(data, 1), size(data, 2));

% if ( nargin < 3 )
%     plotFlag = 0;
% end
% 
% [Nx,Ny,Nc] = size(data);
% 
% %% Sum of squares recon
% rhosos = sqrt(sum(abs(data).^2,3));
% 
% %% Walsh method
% % Block size
% % These have to be odd for the code below
% Bx = kSize(1); By = kSize(2);
% %Bx = 11; By = 11;
% 
% % % Blurred data
% % data2 = zeros(size(data));
% % for n = 1:Nc
% %     data2(:,:,n) = conv2(ones(15,1), ones(15,1), data(:,:,n), 'same');
% % end
% 
% % Usefule shortcut parameters
% bxmin = floor(Bx/2)+1; bymin = floor(By/2)+1;
% bxmax = Nx-floor(Bx/2); bymax = Ny-floor(By/2);
% dx = floor(Bx/2); dy = floor(By/2);
% 
% % Initialize the result
% v1 = zeros(Nx,Ny,Nc);
% s1 = zeros(Nx,Ny);
% u1 = zeros(Nx,Ny);
% v2 = zeros(Nx,Ny,Nc);
% % Loop over the blocks
% % We ignore the edges here for simplicity
% % for x = bxmin:bxmax
% %     for y = bymin:bymax
% for x = 1:Nx
%     for y = 1:Ny
%         
%         indX = x-dx:x+dx;
%         indY = y-dy:y+dy;
%         
%         p = find(indX<1);
%         if ( ~isempty(p) )
%             indX(p) = Nx - indX(p);
%         end
%         
%         p = find(indY<1);
%         if ( ~isempty(p) )
%             indY(p) = Ny - indY(p);
%         end
%         
%         p = find(indX>Nx);
%         if ( ~isempty(p) )
%             indX(p) = -Nx + indX(p);
%         end
%         
%         p = find(indY>Ny);
%         if ( ~isempty(p) )
%             indY(p) = -Ny + indY(p);
%         end        
%         
%         D = reshape(data(indX,indY,:),[Bx*By,Nc]);       
%         [U,S,V] = svd(D);
%         u1(x,y) = mean(U(:,1));
%         v1(x,y,:) = V(:,1)';
%         s1(x,y) = S(1,1);
%         v2(x,y,:) = exp(1i*angle(u1(x,y)))*V(:,1)';
%     end
% end
% %% Coil sensitivities
% mhat  = v1;
% mhat2 = v2;
% 
% csm = v2;
% 
% if ( plotFlag )
%     %% Adaptive combine
%     rhoac  = sum(conj(mhat).*data,3);
%     rhoac2 = sum(conj(mhat2).*data,3);
% 
%     % Display adaptive combine results
%     figure;
%     subplot(2,2,1);
%     imagesc(abs(rhoac)); axis image; axis off; colorbar
%     title('Rho mag');
%     subplot(2,2,2);
%     imagesc(angle(rhoac)); axis image; axis off; colorbar
%     title('Rho phase');
%     subplot(2,2,3);
%     imagesc(abs(rhoac2)); axis image; axis off; colorbar
%     title('Rho2 mag');
%     subplot(2,2,4);
%     imagesc(angle(rhoac2)); axis image; axis off; colorbar
%     title('Rho2 phase');
% end
% 
