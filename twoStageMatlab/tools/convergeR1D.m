function [R,detV,detW] = convergeR1D(X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    sz= size(X);
    n = sz(1); %% size of sample
    d = sz(2); %% dimension
    m = sz(3); %% number of chains
    xbar = mean(X,'all');
    xbarj= reshape(mean(X,1),d,m,[])';
    Bn = cov(xbarj);
    W  = zeros(d);
    for j=1:m
        W = W + cov(X(:,:,j));
    end
    W = (1.0/(double(m)))*W;
    sigmaplus = (double(n-1)/double(n)) * W + Bn;
    V = sigmaplus + Bn/double(m);
%     V = (double(n-1)/double(n)) * W + (1 + 1/double(m)) * Bn;
    lambda1 = max(eig(inv(W)*Bn));
%     nu= 2*V / (V - 1)
    R = (double(n-1)/double(n)) + (double(m+1)/double(m))*lambda1;
%     R = V / W;
    detV = det(V);
    detW = det(W);
end

