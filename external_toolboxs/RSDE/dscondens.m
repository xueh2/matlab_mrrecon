function [Sdata,R_min]=dscondens(data,k)
%
%Multiscale Density Condensation Estimator
%
% Use format: [Sdata,R_min]=dscondens(data,k)
%
% Input:  data:     data set [Nxd]
%         k:        preset number of the samples covered by each condensed point
% Return: Sdata:    selected condensed points
%         R_min:    the radius of the area covered by each condensed point
%
%% Technique reference:
%     P. Mitra, C.A. Murthy and S.K. Pal. "Density based multiscale data condensation",
%     IEEE Transactions on Pattern Analysis and Machine Intelligence, 24:6, 2002.
%
% Copyright Chao He & Mark Girolami
% Last revised on March 11th, 2003
%

Sdata=[];
R_min=[];
[N,d]=size(data);

m=0;
while N>0
    if N<=k
        k=N-1;
    end
    m=m+1;
    fprintf('m=%d  N=%d\n',m,N);
    disc=zeros(N);
    for i=1:N
        dummy=data(i:end,:)-repmat(data(i,:),N-i+1,1);
        d1=sqrt(sum(dummy.^2,2));
        disc(i:N,i)=d1;
        disc(i,i:N)=d1';
    end
    sorted_disc=sort(disc);
    R=sorted_disc(k+1,:);
    R=R';
        
    
    [R_m,I]=min(R);
    R_min=[R_min,R_m];
    Sdata=[Sdata;data(I,:)];
    I1=find(disc(:,I)>2*R_m);
    data=data(I1,:);
    [N,d]=size(data);
end
    

        