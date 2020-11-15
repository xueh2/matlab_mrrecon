%This is a demo to run the Multiscale Density Condensation Estimation.
%
% Usage: run "demods" in Matlab command window.
%
% Copyright Chao He & Mark Girolami
% Last revised on August 22th, 2002


data=load('data.txt');
%Set the parameter of Multiscale Density Condensation
%  which control the reducing rate. 
k=6;

%Multiscale Density 
[Sdata,R_min]=dscondens(data,k);

%Visualisation
%  Points with small red circle are selected condensed points;
%  Large circles indicate the area covered by each selected point.
figure;
plot(data(:,1),data(:,2),'.');
hold on;
plot(Sdata(:,1),Sdata(:,2),'ro');

[r,c]=size(Sdata);
for i=1:r
    minG1=min(data(:,1));
    minG2=min(data(:,2));
    maxG1=max(data(:,1));
    maxG2=max(data(:,2));
    int1=((maxG1-minG1)/50);
    int2=((maxG2-minG2)/50);
    a= minG1:int1:maxG1;
    b= minG2:int2:maxG2;
    [X, Y] = meshgrid(a,b);
    X = X(:);
    Y = Y(:);
    grid = [X Y];
    for j=1:length(grid)
        z(j)=norm(grid(j,:)-Sdata(i,:));
    end
    Z = reshape(z, length(a), length(b));
    contour(a,b, Z, [2*R_min(i),2*R_min(i)]);
end
hold off;
title('Multiscale Density Condensation');