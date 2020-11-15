% This is a demo program to show the selected points obtained by RSDE
%     lie in the centre of the distribution and the shape they form is 
%     somewhat similar to that obtained by Principal Curves.
%
% Usage: run "demoprincipal" in Matlab command window.
%
% Copyright Chao He & Mark Girolami
% Last revised on 5th Jan., 2004

disp('Depending on the speed of your machine this demo might be a bit slow due to the large sample size,');
disp('Please be patient - it is worth the wait...');

H_min=0.05;
H_max=0.5;
K=10;
cross=10;
qp_cg=2;

data=load('principal_curves_data.txt');
Orig_data=data(1:3000,:);
h_opt_parz=0.10;      %obtained by cross-validation
%h_opt_parz=parzen(Orig_data,[H_min,H_max],K,cross)
h_opt_sparse=0.23;    %obtained by minimising ISE between RSDE and Parzen estimates
% h_opt_sparse=rsde_ise(Orig_data,[H_min,H_max],h_opt_parz,qp_cg)
alpha=load('principal_alpha.txt'); %results obtained by the following code.
%[alpha,rsde_est]=rsde(Orig_data,h_opt_sparse,qp_cg);

%visualisation
G=Orig_data; 
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
    P=rsdeprob(Orig_data,grid(i,:),alpha,h_opt_sparse);
    Ps=[Ps;P];
end
Z = reshape(Ps, length(b), length(a)); 

Ps1=[];
for i=1:length(grid)
    P1=parzenprob(Orig_data,grid(i,:),h_opt_parz);
    Ps1=[Ps1;P1];
end
Z1 = reshape(Ps1, length(b), length(a)); 

%Plot the density contour of Parzen window estimation and the original data points
figure
plot_data=Orig_data(1:1000,:); % only 1000 points out of the original dataset are plotted
plot(plot_data(:,1),plot_data(:,2),'b.');
hold on;
contour(a,b,Z1,'r');
hold off;

%Plot the density contour of RSDE and the selected points (encircled by red)
figure
hold on;
I=find(alpha);
for i=1:length(I)
    plot(Orig_data(I(i),1),Orig_data(I(i),2),'bo');
    plot(Orig_data(I(i),1),Orig_data(I(i),2),'b.');
end
contour(a,b,Z,'r');
hold off;
