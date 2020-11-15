function p_x=rsdeprob(Ref_Set,Test_Point,alpha,h)
%
%
%Calculate the probability density of given Test_Point by using Reduced set density estimate obtained by Ref_Set. 
%
%   Use format: p_x=rsdeprob(Ref_Set,Test_Point,alpha,h)
%
%   Input:    Ref_Set:      Training data set;  ([Nxd]: N is the sample amount, d is the dimension);
%             Test_Point:   Testing point;  [1xd]
%             alpha:        Weight vector obtained by RSDE [Nx1]
%             h:            The window width for reduced set
%
%   Return:   p_x [1x1]:    Probability density of testing point fitted by RSDE
%
%   Copyright Mark Girolami & Chao He
%   Last revised on January 23th, 2003
%   Acknowledgement: Thanks to Anna Szymkowiak from Technical University of Denmark 
%                    for vectorising certain parts of the code.


[N,d]=size(Ref_Set);

I=find(alpha); % find non-zero alpha's
Nnz = length(alpha(I)); %NO of non-zero-components

dummy=Ref_Set(I,:)-repmat(Test_Point,Nnz,1);
p_x=sum(alpha(I)./((((2*pi)^(d/2))*h^d)).*exp(-0.5*sum(dummy.^2,2)/(h^2)));
