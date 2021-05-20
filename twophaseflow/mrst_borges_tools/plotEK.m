function plotEK(phi,E,E0)
    m = length(E);
    c1 = sum(log(phi));
    c2 = sum(log(phi.*E));
    c3 = sum(phi);
    sp = (E0*c1 - c2) / (c3 * (2*E0 - 1))
    x = linspace(min(phi),max(phi),1000);
    y = c*(x.^3)./((1-x).^2);
    twoKCplot(phi, x, k, y,...
        '\phi','\kappa','computed','model',c,10);
    pause(3); clf; close all
end

