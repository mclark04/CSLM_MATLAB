function XIT(NAME,N)

if N<0
    error('ERROR: %s %d',NAME,N);
else
    warning('Have a proper shut down here: %s %d',NAME,N);
end
end