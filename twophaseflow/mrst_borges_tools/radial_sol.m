function [p] = radial_sol(x,q,B,mu,rho,rw,h,K,pw)
    p = pw - ((q*mu*B)/(2*pi*K*h))*(log(x/rw));
end

