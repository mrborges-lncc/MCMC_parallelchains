function [out] = prop(theta,mu,s,d)
     out = mvnrnd(theta,s);
%     out = theta + lhsnorm(mu,s,d);
%    out = lhsnorm(mu,s,d);