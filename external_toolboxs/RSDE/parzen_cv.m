function [Av_Vd_Q]=parzen_cv(data,h,cross)

%Parzen Window with cross validation
%[Av_Vd_Q]=parzen_cv(data,h,cross)
%
%   data:         training data set;  [Nxd] (N is the sample amount, d is the dimension);
%   h:            Parzen Window width
%   cross         number of cross validation folds
%
%   Copyright Chao He & Mark Girolami
%   Last revised on August 22th, 2002
%


[N,d]=size(data);

%Cross-Validation
nC=floor(N/cross);
for k=1:cross
    switch k
    case 1
        Vd_data=data(1:nC,:);
        Tr_data=data(nC+1:N,:);
    case cross
        Vd_data=data((cross-1)*nC+1:N,:);
        Tr_data=data(1:(cross-1)*nC,:);
    otherwise
        Vd_data=data((k-1)*nC+1:k*nC,:);
        Tr_data=[data(1:(k-1)*nC,:);data(k*nC+1:N,:)];
    end
    [rV,cV]=size(Vd_data);
    [rT,cT]=size(Tr_data);
    U_data=[Tr_data;Vd_data];      
        
    %Parzen Window
    for i=1:rV
        p_x(i,1)=parzenprob(U_data(1:rT,:),U_data(rT+i,:),h);
    end
    Vd_Q(k)=sum(log(p_x))/rV;
    %fprintf('cross=%d   Vd_Q=%f\n',k,Vd_Q(k));
end

%Calculate the average log-likelihood of training and validation
Av_Vd_Q=sum(Vd_Q)/cross;
