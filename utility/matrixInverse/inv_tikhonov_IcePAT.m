
function psedu_inv_1 = inv_tikhonov_IcePAT(c_0, varargin)
% first order tikhonov regularization using IcePAT strategy
% double dChi = m_pPATFunctor->getGrappaRegularizationWeight(); // 1.0E-4;
% if ( dChi > 1.0E-9 )
% {
%     double dTraceR = 0.0;
%     double dTraceI = 0.0;
%     int n;
%     for ( n=0; n<nRowA; n++ )
%     {
%         dTraceR += (double)mAmAH(n,n).real();
%         dTraceI += (double)mAmAH(n,n).imag();
%     }
%     //   - trace(H) equals nRowA for unity diagonal elements
%     double lambda  = dChi * dTraceR / nRowA;
%     for ( n=0; n<nRowA; n++ )
%     {
%         mAmAH(n,n).set( mAmAH(n,n).real() + (float)lambda, 0.0f );
%     }
%     if (bbTrace) {ICE_OUT("      Regularization: nRowA "<<nRowA<<"; dChi "<<dChi<<"; dTraceR "<<dTraceR<<"; dTraceI "<<dTraceI<<"; lambda "<<lambda<<";");}
% }

psedu_inv_1 = 0;
s = size(c_0);

if nargin == 2
    sigma = varargin{1};
else
    sigma = 0.0001; % 1/10000 of the maximum singular value
end

% if sigma > 0.99999999,
%     'Error! sigma must be a real number (0,1)!'
%     return
% end

if s(1) ~=s(2)
    'Error!, c_0 must be a square matrix'
    return
end
%tic
% if s(1) < 1600
%     [V,D]=eig(single(c_0+(c_0)')/2); %'a', toc, tic,
%     E = abs(diag(D));
%     
%     nRow = size(c_0, 1);
%     lamda  = sigma * sum(real(diag(D))) / nRow;
%     
%     % E = E + sigma*E(end);
%     E = E + lamda;
%     inv_D = diag(1./E);
% else
%     opts.disp=0; opts.tol=0.0001; opts.issym=1; 
%     [V, D] = eigs(double(c_0), 200, 'LM', opts);
%     
%     nRow = size(c_0, 1);
%     lamda  = sigma * sum(real(diag(D))) / nRow;
%     
%     E = abs(diag(D));
%     E = E + lamda;
%     inv_D = diag(1./E);
% end
% psedu_inv_1 = V*inv_D*V' ;

dTraceR = trace(real(c_0));
nRow = size(c_0, 1);
lamda  = sigma * dTraceR / nRow;
for n=1:nRow
    c_0(n,n) = real(c_0(n,n)) + lamda;
end
psedu_inv_1 = inv(c_0);
% psedu_inv_1 = c_0;
