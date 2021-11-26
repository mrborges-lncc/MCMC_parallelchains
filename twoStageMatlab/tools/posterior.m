function [out] = posterior(lk,theta,M)
%     sz  = size(theta,3);
%     mu  = zeros(1,M)
%     s   = eye(M)
    p   = 1.0;
%     for i = 1 : sz
%         p = p * mvnpdf(theta(:,1,i));
%     end
    out = lk * p;
end

