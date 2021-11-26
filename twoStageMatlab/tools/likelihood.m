function [out] = likelihood(cref,csample,sigma,num_datatype)
    out = 0.0;
    for k = 1 : num_datatype
        ref    = cref{k};
        sample = csample{k};
        out    = out - ((norm(ref(:,2:end) - sample(:,2:end))^2)/...
            (norm(ref(:,2:end))^2))/(2*sigma(k));
    end
    out = exp(out);
    return
end