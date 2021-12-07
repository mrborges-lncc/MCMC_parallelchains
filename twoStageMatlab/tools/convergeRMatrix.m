function [R,detV,detW] = convergeRMatrix(X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    sz= size(X);
    n = sz(1);
    d = sz(2);
    m = sz(3);
    xbar = mean(mean(X,1),3);
    xbari= reshape(mean(X,1),d,m,[])';
    Bn = cov(xbari);
    W  = zeros(d,d);
    for j=1:m
        W = W + cov(X(:,:,j));
    end
    W = (1.0/(double(m)))*W;
    V = double((n-1)/n)*W + double(1+1/m)*Bn;
    lambda1 = max(eig(inv(W)*Bn));
    R = double((n-1)/(n)) + double((m+1)/m)*lambda1;
    detV = det(V);
    detW = det(W);
end

