function [out] = erromedio(cref,csample,num_datatype)
    out = zeros(1,num_datatype);
    for k = 1 : num_datatype
        ref    = cref{k};
        sample = csample{k};
        out(1,k) = ((norm(ref(:,2:end) - sample(:,2:end))^2)/...
            (norm(ref(:,2:end))^2));
    end
    return
end