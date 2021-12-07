function [R,V,W] = convergeR(X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
sz= size(X);
n = sz(1);
d = sz(2);
m = sz(3);
xbar = mean(mean(X,1),3);
xbari= reshape(mean(X,1),d,m,[])';
Bn   = zeros(1,d);
for i=1:m
    Bn = Bn + (xbari(i,:)-xbar).^2;
end
Bn  = Bn/double(m-1);
W = 0.0;
for j=1:m
    for t=1:n
        W = W + (X(t,:,j)-xbari(j,:)).^2;
    end
end
W = (1.0/(double(m*(n-1))))*W;
s2= double((n-1)/n)*W+Bn;
V = s2 + Bn/double(m);
R = double((n+3)/(n+1))*V./W;
end

