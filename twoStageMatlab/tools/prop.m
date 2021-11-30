function [out] = prop(method,theta,mu,s,d)
    switch method
        case 'RW'
            out = mvnrnd(theta,s);
        case 'CW'
            beta= s;% * s;
            out = sqrt(1-beta) * theta + sqrt(beta) * lhsnorm(0.0,1.0,d);
        otherwise
            out = lhsnorm(0.0,1.0,d);
    end
%     out = theta + lhsnorm(mu,s,d);
%    out = lhsnorm(mu,s,d);
    return
end