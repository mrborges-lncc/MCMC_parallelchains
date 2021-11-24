function [out] = prop(mu,s,d)
%     out = mvnrnd(mu,s);
    out = lhsnorm(mu,s,d);