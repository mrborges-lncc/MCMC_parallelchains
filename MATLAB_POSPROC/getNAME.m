function [var]=getNAME(name)
n=length(name);
while name(n)~='_',
    n=n-1;
end
var=n;