function [out] = KL(T,theta,KLM,numel)
    out = zeros(numel,1);
    for el = 1:numel
        out(el) = T(el,:)*theta;
    end
end

