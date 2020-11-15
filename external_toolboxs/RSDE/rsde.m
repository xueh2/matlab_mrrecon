function [alpha,rsde_est]=rsde(Ref_Set,h,qp_cg)
%
% Reduced Set Density Estimation (RSDE).
%      
%    Use format: alpha=rsde(Ref_Set,h,qp_cg)
%
%    Input:  Ref_Set [Nxd]:   Training(reference) dateset 
%            h:               Window width for reduced set 
%            qp_cg:           Optimisation method flag
%                             qp_cg=1: Standard Quadratic Optimization
%                             qp_cg=2: Sequential Minimal Optimisation (SMO)
%                             qp_cg=3: Multiplicative update optimisation
%    Return: alpha  [Nx1]:    Weight vector obtained by RSDE
%            rsde_est [1xN]:  Reduced set density estimate for reference dataset
%
%    Copyright Chao He & Mark Girolami
%    Last revised on January 22th, 2003
%    Acknowledgement: Thanks to Anna Szymkowiak from Technical University of Denmark 
%                     for vectorising certain parts of the code.


[N,d]=size(Ref_Set); 
Z=zeros(N, 'single');
%Gaussian pdf kernel is being used
disp('Waiting for calculating kernal functions...');

for i=1:N
  dummy=Ref_Set(i:end,:)-repmat(Ref_Set(i,:),N-i+1,1);
  p_x=-0.5*sum(dummy.^2,2);
  Z(i:N,i)=p_x;
  Z(i,i:N)=p_x';
end

pi2=2*pi;
Q=(1/((pi2^(d/2))*((sqrt(2)*h)^d))).*exp(Z./(2*h^2));
K=(1/((pi2^(d/2))*h^d)).*exp(Z./(h^2));

switch qp_cg
case 1
    %Standard Quadratic Optimization
    A=ones(length(K),length(K));
    b=ones(length(K),1);
    c=-(1/N)*sum(K)';
    alpha=spsolqp(Q,A,b,c);
    %----------------------------
case 2
    %Sequential Minimal Optimisation
    D=sum(K)./N;
    alpha=SMO(Q,D)';
    %---------------------------
case 3
    %Multiplicative update optimisation
    D=sum(K)./N;
    alpha=multupd(Q,D);
    %---------------------------
end
rsde_est=alpha'*K;