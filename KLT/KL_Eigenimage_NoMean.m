% This is to generate the temporal eigenimage 
% [a_eigen,V,D] = KL_Eigenimage(a);
% a: 3-D data set, the first 2 dimensions are spatial, the 3rd is temporal

function [a,V,D, eigImg] = KL_Eigenimage_NoMean(a);

%a = single(a);
s = size(a);

if length(s) ==2, % 2-D case
    if s(2) > s(1), a = a'; end % first dimension is always larger
    s = size(a);
    t_mean = mean(a,1);
    a_mean = ones(s(1),1)*t_mean;
    a = a - a_mean; %remove mean
    A = (a'*a)/(s(1)); %clear b;
    [V,D] = (eig((A+A')/2));
    %[V,D] = eig((A+A')/2, diag(ones(1,s(2))), 'chol'); % Diagonalize the covariance matrix
    eigImg = a*V;
    a = (a + a_mean)*V;
elseif length(s)==3, % 3-D case
    % b = single(zeros(N,s(1)*s(2)));
    a = (reshape(a, [s(1)*s(2),s(3)] )); %clear a
    t_mean = mean(a,1);
    a_mean = ones(s(1)*s(2),1)*t_mean;
    a = a - a_mean; %remove mean
    A = (a'*a)/(s(1)*s(2)); %clear b;
    [V,D] = (eig((A+A')/2)); % Diagonalize the covariance matrix
    %[V,D] = eig((A+A')/2, diag(ones(1,s(3))), 'chol'); 
    V = (V);
    D = (D);
    eigImg = (reshape( a*V , [s(1),s(2),s(3)] ) ) ;;
    a = (reshape( (a+a_mean)*V , [s(1),s(2),s(3)] ) ) ; clear b
else
    'Input array is not a 2/3-D array!'
end


% % Old one
% if length(s) ~=3,
%     'Input array is not a 3-D array!'
% else
%     %b = single(zeros(N,s(1)*s(2)));
%     a = (reshape(a, [s(1)*s(2),s(3)] )); %clear a
%     t_mean = mean(a,1);
%     a_mean = ones(s(1)*s(2),1)*t_mean;
%     a = a - a_mean; %remove mean
%     A = (a'*a)/(s(1)*s(2)); %clear b;
%     [V,D] = (eig((A+A')/2)); % Diagonalize the covariance matrix
%     V = (V);
%     D = (D);
%     a = (reshape( (a+a_mean)*V , [s(1),s(2),s(3)] ) ) ; clear b
% 
% end










