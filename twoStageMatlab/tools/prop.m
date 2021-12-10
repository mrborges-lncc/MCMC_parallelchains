function [out] = prop(method,theta,chain,nk,s,d,NC,freqj,iter)
    if mod(iter,freqj) == 0, s = 2 * s; end
    if iter < 100, s = 10 * s; end
    switch method
        case 'RW'
            out = mvnrnd(theta(:,chain,nk),s);
        case 'CW'
            out = sqrt(1-s^2) * theta(:,chain,nk) + s * lhsnorm(0.0,1.0,d);
        case 'DE'
            [r] = chooser(chain,NC);
            out = theta(:,chain,nk) + ...
                s*(theta(:,r(1),nk) - theta(:,r(2),nk)) + lhsnorm(0.0,1.0e-03,d);
        otherwise
            out = lhsnorm(0.0,1.0,d);
    end
    return
end