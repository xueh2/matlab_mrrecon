function [v,d]=eig_power(R);
% function [v,d]=eig_power(R);
%
% vectorized method for calculating the dominant eigenvector based on
% power method. Input, R, is an image of sample correlation matrices
% where: R(y,x,:,:) are sample correlation matrices (ncoil x ncoil) for each pixel
%
% v is the dominant eigenvector
% d is the dominant (maximum) eigenvalue

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

rows=size(R,1);cols=size(R,2);ncoils=size(R,3);
N_iterations=2;
v=ones(rows,cols,ncoils); % initialize e.v.

d=zeros(rows,cols);
for i=1:N_iterations
    v=squeeze(sum(R.*repmat(v,[1 1 1 ncoils]),3));
	d=rss(v);
	v=v./repmat(d,[1 1 ncoils]);
end

% % try to use the eigen image for phase removal
[eigIm,eigV,eigD] = KL_Eigenimage(v);
% 
% % used mean 
% vUsed = mean(v, 3);
% 
vUsed = v(:,:,1);
% 
% % vUsed = sum(eigIm, 3);
% 
% vUsed = eigIm(:,:,end);


vUsed(find(abs(vUsed)==0)) = eps;
vUsed = vUsed ./ (abs(vUsed));
 
v = v.*repmat(conj(vUsed),[1 1 ncoils]);

v=conj(v);

return
