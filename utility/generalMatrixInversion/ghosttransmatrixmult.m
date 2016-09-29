function [st,gmat] = ghosttransmatrixmult(s,tformdata,transp_flag)
% ST = GHOSTTRANSMATRIXMULT(S,TFORMDATA,TRANSP_FLAG)
%
% ST is the image S ghosted as described by tformdata
% ST is passed as a N by array, where N = n1*n2*...nd for a d-dim
% image.
%
%    a structure with fields:
% 
%      TFORMDATA.TFORMS: array of Matlab image transforms
%      TFORMDATA.SHOTS:  array of positions in k-space (same length)
%      TFORMDATA.SIZ:    size of input image
%      
%      TRANSP_FLAG is optional, used by iterative algorithms such as 'LSQR'
%      if set to 'transp': returns the transpose operation
%
%    this means:
%           
%       -at time t, s is transformed according to tforms{t}
%       -the transformed image is Fourier transformed (FT)
%       -the positions of the FT corresponding to shots{t} are kept
%       -st is the IFT (inverse FT) of the mixed Fourier transforms
%       
% This corresponds to multiplication of s(:) with a 'ghosting' matrix.
% As the aim is to invert it, the function also allows the
% TRANSP_FLAG to be set to 'transp', as required by  
% Matlab CG type iterative solvers 
%   
% SHOTS{T} must be linear indices, this allows any dimension.
%
%
% [ST,GMAT] = GHOSTTRANSMATRIXMULT(S,TFORMDATA,TRANSP_FLAG)
% returns a sparse matrix  such that GMAT*S(:) '=' ST
% or GMAT'*S(:) if 'transp' is set.
% Sparsity is controlled by TFORMDATA.SPTHRESH 
% whose default is (max magnitude of s)/10000. 
% 
%
% WARNING: 1) The building of the sparse matrix is very slow, as it uses
% the 'pedestrian' method of applying the operation to images of a point.
% The functions GHOSTTRANSMATRIX try to be more efficient, but for the moment
% there is a problem (the output appears to have different phases).
% This matrix version can be useful if the motions at each shot are known, say e.g.
% from some navigator/tagged registrations, but the inverse transform is not, a common problem
% with non-rigid registrations. 
%          2) The code uses the Image Processing Toolbox.
% 
% Computation times are:
% < ~ 0.5s for 96x96 image output, thus ~96*96*0.5s for sparse mat. 
%
% REFERENCES
% 
% @Article{batchelor05:_matrix_descr_of_gener_motion,
%   author = 	 {P. G. Batchelor and  D. Atkinson and P. Irarrazaval and D. L. G. Hill and J. Hajnal and D. Larkman },
%   title = 	 {Matrix description of general motion correction applied to multishot images},
%   journal = 	 {Magn. Res. Med.},
%   year = 	 2005,
%   volume = 	 54,
%   number = 	 5,
%   pages = 	 {1273--1280}
% }
% 
% 
%
% SEE ALSO:
% LSQR, PCG, GHOSTTRANSMATRIX, IPT
%
%------------------------------------------------------- 
% <EXAMPLE (Need image processing toolbox)>
% s0 = phantom;
%
% % load any  image, call it s0.
% n1 = size(s0,1); n2 = size(s0,2);
% [X,Y] = meshgrid(1:n2,1:n1);
% X0 = n2/2; Y0 = n1/2;
% tr0 = [n2,n1];
% r = min(X0,Y0)/3;
% % random affine transforms in 4 shots
% nshots = 4;
% S = zeros(size(s0));
% A0 = eye(3);
% for t = 1:nshots
%   A = A0 + 0.1*[randn(3,2),[0;0;1]]; A(3,3) = 1; % random transform
%   tforms{t} = maketform('affine',A);
%   shots{t} = sub2ind(size(s0),Y(t:nshots:end),X(t:nshots:end));
%
%   st = imtransform(s0,tforms{t},'cubic','xdata',[1,n2],'ydata',[1,n1]);
%
%   St = fft2(st);
%   S(shots{t}) = St(shots{t});
% end
% s = ifft2(S);
% % INPUT FOR ghosttransmatrixmult
% tfd.tforms = tforms;
% tfd.shots  = shots;
% tfd.interpolation = 'cubic';
% tfd.siz = size(s0);
% maxit = 15;
% sh = lsqr(@ghosttransmatrixmult,s(:),[],maxit,[],[],[],tfd);
% sh = reshape(sh,size(s0));
%
% </EXAMPLE>
%-------------------------------------------------------
%
% ACKNOWLEDGMENTS
% K Pruessman.
% EPSRC GR/S30184 and AF/001381
% Chilean FONDECYT 1030570.
% 
% (c) philip.batchelor@ucl.ac.uk

tforms = tformdata.tforms;
shots  = tformdata.shots;
if isfield(tformdata,'interpolation')
  R      = tformdata.interpolation;
else
  warning('using linear interpolation');  
  R      = 'linear';
end

siz    = tformdata.siz;
s = reshape(s,siz);
nshots = length(tforms);

if ndims(s) == 2
    
  xd = [1,size(s,2)];
  yd = [1,size(s,1)];
  %hs = waitbar(0,'Shots processed...');
  if (nargin > 2) & strcmp(transp_flag,'transp')
    st = zeros(size(s));
    for t = 1:nshots
      %waitbar(t/nshots,hs);
      S = zeros(size(s));
      St = fft2(s);
      S(shots{t}) = St(shots{t});
      sa = ifft2(S);
      st = st + imtransform(sa,fliptform(tforms{t}),R,'xdata',xd,'ydata',yd);
    end
  else
    S = zeros(size(s));
    for t = 1:nshots
      %waitbar(t/nshots,hs);
      St = fft2(imtransform(s,tforms{t},R,'xdata',xd,'ydata',yd));
      S(shots{t}) = St(shots{t});
    end
    st = ifft2(S);
  end
  %close(hs);

else
  if (nargin > 2) & strcmp(transp_flag,'transp')
    st = zeros(size(s));
    for t = 1:nshots
      S = zeros(size(s));
      St = fftn(s);
      S(shots{t}) = St(shots{t});
      sa = ifftn(S);
      st = st + tformarray(sa,fliptform(tforms{t}),R,[1,2,3],[1,2,3],size(s),[],0);
    end
  else
    S = zeros(size(s));
    for t = 1:nshots
      St = fftn(tformarray(s,tforms{t},R,[1,2,3],[1,2,3],size(s),[],0));
      S(shots{t}) = St(shots{t});
    end
    st = ifftn(S);
  end
end
st = st(:);
%=======================================================
% fftw('planner','exhaustive');
% St = fft2(s);
% str = fftw('wisdom');


if nargout == 2

  if ~isfield(tformdata,'spthresh')
    spthresh = 1.0e-4*max(abs(s(:)));
  else
    spthresh = tformdata.spthresh*max(abs(s(:)));
  end

  fnz    = 15; % preallocate 15 non-zero element/pixel
	       % but will resize if nec
	       % could estimate from interpolation?

   [n1,n2] = size(s); N = numel(s);
   O     = zeros(n1,n2);
   I2 = zeros(fnz*N,1); J2 = I2; V2 = I2;
   s1 = zeros(n1,n2);
   l = 1;
   h = waitbar(0,'Pixels processed...');

   for k = 1:N
     waitbar(k/N,h);
     s1(k) = 1;
     if k > 1; s1(k-1) = 0; end; 
     S = O;
      for t = 1:nshots
       St = fft2(imtransform(s1,tforms{t},'cubic','xdata',xd,'ydata',yd));
       S(shots{t}) = St(shots{t});
      end
     s2 = ifft2(S);  

     if isfield(tformdata,'disp_psf')
        if tformdata.disp_psf
         imshow(mat2gray(abs(s2)),[],'InitialMagnification','fit'); drawnow;
        end
     end
     [i2,j2] = find(abs(s2) > spthresh);
     ii = sub2ind(size(s1),i2(:),j2(:));
     nz = length(ii);
     jj = repmat(k,nz,1);
     I2(l : l + nz - 1) = ii(:);
     J2(l : l + nz - 1) = jj(:);
     V2(l : l + nz - 1) = s2(ii(:)); %v2(:);
     l = l + nz;
   end
   % free unused memory
   if l+nz+1 < length(I2) %+nz??
     I2(l+nz:end) = []; J2(l+nz:end) = []; V2(l+nz:end) = [];
   end
   
   if (nargin > 2) & strcmp(transp_flag,'transp')
     gmat = sparse(J2,I2,V2,N,N);
   else	 
     gmat = sparse(I2,J2,V2,N,N);	   
   end
   close(h);
end


