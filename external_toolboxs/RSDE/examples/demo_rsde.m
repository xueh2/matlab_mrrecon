% This is a demo program to compare the density estimation results respectively 
%    obtained by the Parzen Window estimator and the Reduce Set Density Estimator
%    (RSDE) for 2-dimensional reference data. The 2-D visualisation of the results
%    is also given.
% 
% Usage: run "demo_rsde" in Matlab command window.
%
% Change the flag "qp_cg" to select the different optimisation methods.
%     qp_cg=1:  Standard quadratic optimization
%     qp_cg=2:  Sequential minimal optimisation(SMO) - suggested to use, fastest
%     qp_cg=3:  Multiplicative update optimisation
%
% Copyright Chao He & Mark Girolami
% Last revised on 5th Jan., 2004

qp_cg=2;


%Parameters for using function "parzen", see help in parzen.m for their meanings
H_min=0.1;
H_max=1;
K=20;
cross=10;

Ref_data=load('data.txt');     

%-----------------------------------------------------------------------------

disp('Parzen window density estimation...\n');

h_opt_parz=parzen(Ref_data,[H_min,H_max],K,cross);
%result: h_opt_parz=0.43;
h_opt_parz=0.43;

fprintf('\t h_opt_parz=%f\n\n',h_opt_parz);

%-----------------------------------------------------------------------------

disp('Reduced set density estimation...\n');

%h_opt_sparse=rsde_ise(Ref_data,[H_min,H_max],h_opt_parz,qp_cg); 
%result: h_opt_sparse=0.73;
h_opt_sparse=0.73;

fprintf('\t h_opt_sparse=%f\n\n',h_opt_sparse);
[alpha,rsde_est]=rsde(Ref_data,h_opt_sparse,qp_cg);

%------------------------------------------------------------------------------

%Visualisation
G=Ref_data; 
miG1=min(G(:,1));miG2=min(G(:,2));
maG1=max(G(:,1));maG2=max(G(:,2));
int1=((maG1-miG1)/50);
int2=((maG2-miG2)/50);
a= miG1:int1:maG1;
b= miG2:int2:maG2;

[X, Y] = meshgrid(a,b);
X = X(:);
Y = Y(:);
grid = [X Y];
   
Ps=[];
for i=1:length(grid)
    P=rsdeprob(Ref_data,grid(i,:),alpha,h_opt_sparse);
    Ps=[Ps;P];
end      
Z = reshape(Ps, length(b), length(a)); 
   
Pg=[];
for i=1:length(grid)
  	P=parzenprob(Ref_data,grid(i,:),h_opt_parz);
  	Pg=[Pg;P];
end      
ZZ = reshape(Pg, length(b), length(a)); 

%Plot the density contour obtained by Parzen density estimation and the reference data points
figure
contour(a,b,ZZ,'r');
hold on;
plot(Ref_data(:,1),Ref_data(:,2),'b.');
title('Parzen Window Density Estimation')
hold off;

%Plot the density contour obtained by RSDE and the selected points
figure
I=find(alpha); % if using standard quadratic optimisation, i.e. qp_cg=1, threshold for alpha will be needed then,
               % e.g. I=find(alpha>0.00001).
contour(a,b,Z,'r');
hold on;
plot(Ref_data(I,1),Ref_data(I,2),'bo');
plot(Ref_data(I,1),Ref_data(I,2),'b.');
title('Reduced Set Density Estimation')
hold off;

h_opt_parz
h_opt_sparse