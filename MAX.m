function [var]=MAX(var,maxVar)

if nargout<1
    error('Need to pass out the min data');
end

tmp=var;
var(var<maxVar)=maxVar;