function h_best=rsde_ise(Ref_Set,H,h_parz,qp_cg)
%
%Minimising the ISE between the RSDE and the existing Parzen estimate to give the optimised window width value for RSDE.
%
%   Use format: h_best=rsde_ise(Ref_Set,H,h_parz,qp_cg)
%
%   Input:    Ref_Set [N x d]:   Reference data set
%             H       [1 x 2]:   [H_min,H_max] denotes the searching range of the RSDE window width
%             h_parz:            Parzen window width
%             qp_cg:             Optimisation method flag
%                                qp_cg=1: Standard Quadratic Optimization
%                                qp_cg=2: Sequential Minimal Optimisation (SMO)
%                                qp_cg=3: Multiplicative update optimisation
%
%   Return:   h_best:            Optimised RSDE window width
%
%   Copyright Mark Girolami & Chao He
%   Last revised on January 5th, 2004
%   Acknowledgement: Thanks to Anna Szymkowiak from Technical University of Denmark 
%                    for vectorising certain parts of the code.

[N,d]=size(Ref_Set);
KL=[];
INTERVAL=(H(2)-H(1))/10;
K=zeros(N);
%Gaussian pdf kernel is being used

for i=1:N
  dummy=Ref_Set(i:end,:)-repmat(Ref_Set(i,:),N-i+1,1);
  p_x=((1/(((2*pi)^(d/2))*h_parz^d))*exp(-0.5*sum(dummy.^2,2)/(h_parz^2)));
  K(i:N,i)=p_x;
  K(i,i:N)=p_x';
end

parz_est=mean(K);

for h=H(1):INTERVAL:H(2)
   [alpha,rsde_est]=rsde(Ref_Set,h,qp_cg);
   kl=mean((parz_est-rsde_est).^2./parz_est);  %updated on 5th Jan., 2004 by discovery of Alexander Ihler
   fprintf('\t KL=%6.4f h=%6.4f\n',kl, h);
   KL=[KL;kl h]; 
end

[a,b]=min(KL(:,1));
h_best=KL(b,2);
%fprintf('\t KL_min=%6.4f  h_best=%6.4f\n\n', a, h_best);
%plot(KL(:,2),KL(:,1))