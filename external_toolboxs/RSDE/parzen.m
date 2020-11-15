function h_best=parzen(data,H,K,cross)


%Fit a Parzen Window to the training data
%by using a validate set to select the best parameters.

%   Use format:  h_best=parzen(data,H,K,cross)
%
%   Input:  data:           training data set;  NxP (N is the sample amount, P is the dimension);
%           H:              [H_min,H_max] denotes the searching range of the Parzen window width;
%           K:              total steps intended to search in H;
%           cross:          number of cross validation folds;
%   return: h_best:         Parameter of Parzen window density estimator.
%
%   Copyright Chao He & Mark Girolami
%   Last revised on August 22th, 2002
%

%fprintf('cross validation...\n\n');
h_delta=(H(2)-H(1))/(K-1);
for i=1:K
    fprintf('\t i=%d\t',i);
    h(i)=H(1)+h_delta*(i-1);
    Av_Vd_Q(i)=parzen_cv(data,h(i),cross);
    fprintf('h=%f, Av_Vd_Q=%f\n',h(i),Av_Vd_Q(i));
end
%figure(1);
%i=1:K;
%plot(h(i),Av_Vd_Q(i),'r');
%title('cross validation');
[Q_max,Q_i]=max(Av_Vd_Q);
h_best=h(Q_i);
%fprintf('h=%d\n\n',h_best);