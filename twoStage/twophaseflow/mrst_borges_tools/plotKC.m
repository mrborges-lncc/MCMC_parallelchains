function plotKC(phi,k,ckc)
    m = length(phi);
    c = sum(k.*(phi.^3)./((1-phi).^2)) / sum((phi.^6)./((1-phi).^4));
    x = linspace(min(phi),max(phi),1000);
    y = c*(x.^3)./((1-x).^2);
    twoKCplot(phi, x, k, y,...
        '\phi','\kappa','computed','model',c,10);
end

