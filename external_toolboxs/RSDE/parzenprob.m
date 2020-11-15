function p_x=parzenprob(Ref_Set,Test_Point,h)
%
%
%Calculate the probability density of given Test_Point by using Parzen window density estimate obtained by Ref_Set. 
%
%   Use format: p_x=parzenprob(Ref_Set,Test_Point,h)
%
%   Input:    Ref_Set:      Training data set;  ([Nxd]: N is the sample amount, d is the dimension);
%             Test_Point:   Testing point;  [1xd]
%             h:            The window width for Parzen Window
%
%   Return:   p_x [1x1]:    Probability density of testing point fitted by Parzen window
%
%   Copyright Mark Girolami & Chao He
%   Last revised on January 22th, 2003
%   Acknowledgement: Thanks to Anna Szymkowiak from Technical University of Denmark 
%                    for vectorising certain parts of the code.

[N,d]=size(Ref_Set);
 
dummy=Ref_Set-repmat(Test_Point,N,1);
p_x=sum((1/(N*((2*pi)^(d/2))*h^d))*exp(-0.5*sum(dummy.^2,2)/(h^2)));
