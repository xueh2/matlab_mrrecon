function y = groupSoftThresholding3(x,mu)
	%group softresholding of variable x with treshold mu
	%the groups are formed along the 3rd dimension
    num=size(x,3);
	n=sqrt(sum(x.*conj(x),3));
	w=n-mu;
    w=(abs(w)+w)/2;
    i=(w~=0);
	w(i)=w(i)./n(i);
    y=x.*repmat(w,[1 1 num]);
	
end
