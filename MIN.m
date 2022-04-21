function [var]=MIN(var,minVar)

if nargout<1
    error('Need to pass out the min data');
end

tmp=var;
var(var>minVar)=minVar;