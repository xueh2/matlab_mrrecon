function out=clamp(in, value);
% function out=clamp(in, value);
%
% function will clamp the input to "value" for values > "value"


out=in;

if nargin<2
    return
elseif isempty(value)
    return
elseif length(value)==1 % assume value is maxvalue for clamp
    out(in > value)=value;
elseif length(value)==2 % assume value is vector of min and max values for clamp
    out(in > value(2)) = value(2);
    out(in < value(1)) = value(1);
end
    
    
    