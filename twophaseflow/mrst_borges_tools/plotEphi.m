function plotEphi(phi,E,E0)
    m = length(E);
    syms f a
    f(a) = E0 * exp(-a .* phi);
    lnf(a) = log(f);
    syms err(a) derr
    err(a) = sum(log(E) - lnf)^2;
    derr =  diff(err,a) == 0;
    sp= double(solve(derr,a));
    x = linspace(min(phi),max(phi),1000);
    y = E0 * exp(-sp * x);
    twoEphiplot(phi, x, E, y,...
        '\phi','\mathsf{E}','computed','model',E0, -sp, 11);
end

